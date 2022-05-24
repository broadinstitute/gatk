package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.tribble.Feature;
import htsjdk.variant.vcf.VCFHeader;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.filters.CountingReadFilter;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.tools.sv.SVFeaturesHeader;

import java.util.*;

/**
 * A MultiFeatureWalker is a tool that presents one {@link Feature} at a time in sorted order from
 * multiple sources of Features.  The input files for each feature must be sorted by locus.
 *
 * To use this walker you need only implement the abstract
 * {@link #apply(F, Object, ReadsContext, ReferenceContext)} method in a class that declares
 * a collection of FeatureInputs as an argument.
 * You may optionally implement {@link #onTraversalStart()}, {@link #onTraversalSuccess()},
 * and/or {@link #closeTool()}.
 */
public abstract class MultiFeatureWalker<F extends Feature> extends WalkerBase {

    private SAMSequenceDictionary dictionary;
    private final Set<String> samples = new TreeSet<>();

    @Override
    public boolean requiresFeatures(){
        return true;
    }

    @Override
    public String getProgressMeterRecordLabel() { return "features"; }

    @Override
    void initializeFeatures() {
        features = new FeatureManager(this, FeatureDataSource.DEFAULT_QUERY_LOOKAHEAD_BASES,
                            cloudPrefetchBuffer, cloudIndexPrefetchBuffer, getGenomicsDBOptions());
    }

    /**
     * Operations performed just prior to the start of traversal.
     */
    @Override
    public final void onStartup() {
        super.onStartup();
        setDictionaryAndSamples();
    }

    /**
     * {@inheritDoc}
     *
     * Implementation of Feature-based traversal.
     *
     * NOTE: You should only override {@link #traverse()} if you are writing a new walker base class
     * in the engine package that extends this class. It is not meant to be overridden by tools
     * outside the engine package.
     */
    @Override
    public void traverse() {
        CountingReadFilter readFilter = makeReadFilter();
        final MergingIterator<F> iterator =
                new MergingIterator<>(dictionary, features, userIntervals);
        while ( iterator.hasNext() ) {
            final PQEntry<F> entry = iterator.next();
            final F feature = entry.getFeature();
            final SimpleInterval featureInterval = new SimpleInterval(feature);
            apply(feature,
                    entry.getHeader(),
                    new ReadsContext(reads, featureInterval, readFilter),
                    new ReferenceContext(reference, featureInterval));
            progressMeter.update(feature);
        }
    }

    /**
     * Process an individual feature.
     * In general, subclasses should simply stream their output from apply(), and maintain as little
     * internal state as possible.
     *
     * @param feature Current Feature being processed.
     * @param header Header object for the source from which the feature was drawn (may be null)
     * @param readsContext An object that allows querying for the reads the overlap the feature
     * @param referenceContext An object that allows querying for the reference sequence associated with the feature
     */
    public abstract void apply( final F feature,
                                final Object header,
                                final ReadsContext readsContext,
                                final ReferenceContext referenceContext );

    /**
     * Get the dictionary we settled on
     */
    public SAMSequenceDictionary getDictionary() { return dictionary; }

    /**
     * Get the list of sample names we accumulated
     */
    public Set<String> getSampleNames() { return Collections.unmodifiableSet(samples); }

    /**
     * Choose the most comprehensive dictionary available (see betterDictionary method below),
     * and concatenate the sample names available from each feature input.
     * Each feature input may have its own dictionary, and the user can specify an additional master
     * dictionary, reference dictionary, and reads dictionary.  This method makes certain that all
     * of these dictionaries are consistent with regard to contig name and order.  It's OK if one
     * dictionary is a subset of another:  we'll choose the most comprehensive dictionary.
     * (Can't use the getBestAvailableSequenceDictionary method -- it throws if there are multiple
     * dictionaries available.)
     */
    private void setDictionaryAndSamples() {
        DictSource dictSource = new DictSource(getMasterSequenceDictionary(),
                                        StandardArgumentDefinitions.SEQUENCE_DICTIONARY_NAME);
        if ( hasReference() ) {
            final DictSource refDictSource = new DictSource(reference.getSequenceDictionary(),
                                                  StandardArgumentDefinitions.REFERENCE_LONG_NAME);
            dictSource = betterDictionary(refDictSource, dictSource);
        }
        if ( hasReads() ) {
            final DictSource readsDictSource = new DictSource(reads.getSequenceDictionary(), "read-source");
            dictSource = betterDictionary(readsDictSource, dictSource);
        }
        for ( final FeatureInput<? extends Feature> input : features.getAllInputs() ) {
            final Object header = features.getHeader(input);
            if ( header instanceof SVFeaturesHeader ) {
                final SVFeaturesHeader svFeaturesHeader = (SVFeaturesHeader)header;
                final DictSource featureDictSource = new DictSource(svFeaturesHeader.getDictionary(),
                                                                    input.getName());
                dictSource = betterDictionary(featureDictSource, dictSource);
                final List<String> sampleNames = svFeaturesHeader.getSampleNames();
                if ( sampleNames != null ) {
                    samples.addAll(svFeaturesHeader.getSampleNames());
                }
            } else if (header instanceof VCFHeader ) {
                final VCFHeader vcfHeader = (VCFHeader)header;
                final DictSource featureDictSource = new DictSource(vcfHeader.getSequenceDictionary(),
                                                                        input.getName());
                dictSource = betterDictionary(featureDictSource, dictSource);
                samples.addAll(vcfHeader.getSampleNamesInOrder());
            }
        }
        if ( dictSource.getDictionary() == null ) {
            throw new UserException("No dictionary found.  Provide one as --" +
                    StandardArgumentDefinitions.SEQUENCE_DICTIONARY_NAME + " or --" +
                    StandardArgumentDefinitions.REFERENCE_LONG_NAME + ".");
        }
        dictionary = dictSource.getDictionary();
    }

    /**
     * Makes sure that the two dictionaries are consistent with regard to contig names and order.
     * Returns the more comprehensive (larger) dictionary if they're consistent.
     */
    private static DictSource betterDictionary( final DictSource newDict,
                                                final DictSource curDict ) {
        if ( curDict.getDictionary() == null ) return newDict;
        if ( newDict.getDictionary() == null ) return curDict;
        final DictSource smallDict;
        final DictSource largeDict;
        if ( newDict.getDictionary().size() <= curDict.getDictionary().size() ) {
            smallDict = newDict;
            largeDict = curDict;
        } else {
            smallDict = curDict;
            largeDict = newDict;
        }
        int lastIdx = -1;
        final SAMSequenceDictionary largeDictionary = largeDict.getDictionary();
        for ( final SAMSequenceRecord rec : smallDict.getDictionary().getSequences() ) {
            final int newIdx = largeDictionary.getSequenceIndex(rec.getContig());
            if ( newIdx == -1 ) {
                throw new UserException("Contig " + rec.getContig() + " in the dictionary read from " +
                        smallDict.getSource() + " does not appear in the larger dictionary read from " +
                        largeDict.getSource());
            }
            if ( newIdx <= lastIdx ) {
                final String prevContig = largeDictionary.getSequence(lastIdx).getContig();
                throw new UserException("Contigs out of order: Contig " + rec.getContig() +
                        " comes before contig " + prevContig + " in the dictionary read from " +
                        largeDict.getSource() + ", but follows it in the dictionary read from " +
                        smallDict.getSource());
            }
            lastIdx = newIdx;
        }
        return largeDict;
    }

    public static final class DictSource {
        private final SAMSequenceDictionary dictionary;
        private final String source;

        public DictSource( final SAMSequenceDictionary dictionary, final String source ) {
            this.dictionary = dictionary;
            this.source = source;
        }

        public SAMSequenceDictionary getDictionary() { return dictionary; }
        public String getSource() { return source; }
    }

    public static final class MergingIterator<F extends Feature> implements Iterator<PQEntry<F>> {
        final SAMSequenceDictionary dictionary;
        final PriorityQueue<PQEntry<F>> priorityQueue;

        @SuppressWarnings("unchecked")
        public MergingIterator( final SAMSequenceDictionary dictionary,
                                final FeatureManager featureManager,
                                final List<SimpleInterval> intervals ) {
            final Set<FeatureInput<? extends Feature>> inputs = featureManager.getAllInputs();
            this.dictionary = dictionary;
            this.priorityQueue = new PriorityQueue<>(inputs.size());
            for ( final FeatureInput<? extends Feature> input : inputs ) {
                final Iterator<F> iterator =
                        (Iterator<F>)featureManager.getFeatureIterator(input, intervals);
                final Object header = featureManager.getHeader(input);
                addEntry(new PQContext<>(iterator, dictionary, header));
            }
        }

        @Override
        public boolean hasNext() {
            return !priorityQueue.isEmpty();
        }

        @Override
        public PQEntry<F> next() {
            final PQEntry<F> entry = priorityQueue.poll();
            if ( entry == null ) {
                throw new NoSuchElementException("iterator is exhausted");
            }
            final PQEntry<F> newEntry = addEntry(entry.getContext());
            if ( newEntry != null ) {
                if ( newEntry.compareTo(entry) < 0 ) {
                    final Feature feature = newEntry.getFeature();
                    throw new UserException("inputs are not sorted at " +
                                            feature.getContig() + ":" + feature.getStart());
                }
            }
            return entry;
        }

        private PQEntry<F> addEntry( final PQContext<F> context ) {
            final Iterator<F> iterator = context.getIterator();
            if ( !iterator.hasNext() ) {
                return null;
            }
            final PQEntry<F> entry  = new PQEntry<>(context, iterator.next());
            priorityQueue.add(entry);
            return entry;
        }
    }

    public static final class PQContext<F extends Feature> {
        private final Iterator<F> iterator;
        private final SAMSequenceDictionary dictionary;
        private final Object header;

        public PQContext( final Iterator<F> iterator,
                          final SAMSequenceDictionary dictionary,
                          final Object header ) {
            this.iterator = iterator;
            this.dictionary = dictionary;
            this.header = header;
        }

        public Iterator<F> getIterator() { return iterator; }
        public SAMSequenceDictionary getDictionary() { return dictionary; }
        public Object getHeader() { return header; }
    }

    public static final class PQEntry<F extends Feature> implements Comparable<PQEntry<F>> {
        private final PQContext<F> context;
        private final int contigIndex;
        private final F feature;

        public PQEntry( final PQContext<F> context, final F feature ) {
            this.context = context;
            this.feature = feature;
            this.contigIndex = context.getDictionary().getSequenceIndex(feature.getContig());
        }

        public PQContext<F> getContext() { return context; }
        public F getFeature() { return feature; }
        public Object getHeader() { return context.getHeader(); }

        @Override
        public int compareTo( PQEntry<F> entry ) {
            int result = Integer.compare(contigIndex, entry.contigIndex);
            if ( result == 0 ) {
                result = Integer.compare(feature.getStart(), entry.feature.getStart());
                if ( result == 0 ) {
                    result = Integer.compare(feature.getEnd(), entry.feature.getEnd());
                }
            }
            return result;
        }
    }
}
