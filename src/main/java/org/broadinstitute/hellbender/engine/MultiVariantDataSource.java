package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.MergingIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextComparator;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SequenceDictionaryUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import org.broadinstitute.hellbender.utils.variant.VcfUtils;

import java.nio.file.Path;
import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;

/**
 * MultiVariantDataSource aggregates multiple FeatureDataSources of variants, and enables traversals and queries
 * over those sources through a single interface.
 *
 * Two basic operations are available on this data source:
 *
 * -Iteration over all Variants in the data sources, optionally restricted to Variants overlapping
 *  a set of intervals if intervals are provided via {@link #setIntervalsForTraversal(List)}. Traversal
 *  by a set of intervals requires the files to have been indexed using the bundled tool IndexFeatureFile.
 *  The set of intervals provided MUST be non-overlapping and sorted in increasing order of start position.
 *
 * -Targeted queries by one interval at a time. This also requires the files to have been indexed using
 *  the bundled tool IndexFeatureFile. Targeted queries by one interval at a time are unaffected by
 *  any intervals for full traversal set via {@link #setIntervalsForTraversal(List)}.
 */
public final class MultiVariantDataSource implements GATKDataSource<VariantContext>, AutoCloseable {
    private static final Logger logger = LogManager.getLogger(MultiVariantDataSource.class);

    /**
     * List of FeatureDataSource objects aggregated by this MultiVariantDataSource
     */
    final private List<FeatureDataSource<VariantContext>> featureDataSources = new ArrayList<FeatureDataSource<VariantContext>>();

    /**
     * Merged VCF header used for this (aggregate) source, derived from the individual sources.
     */
    final private VCFHeader mergedHeader;

    /**
     * Iterator representing an open traversal over this data source initiated via a call to {@link #iterator}
     * (null if there is no open traversal). We need this to ensure that each iterator is properly closed,
     * and to enforce the constraint (required by Tribble) that we never have more than one iterator open
     * over our feature reader.
     */
    private CloseableIterator<VariantContext> currentIterator;
    private SortedSet<String> mergedSamples;

    /**
     * Creates a MultiVariantDataSource backed by the provided FeatureInputs. We will look ahead the specified number of bases
     * during queries that produce cache misses.
     *
     * @param featureInputs List of FeatureInput<VariantContext>> specifying sources of VariantContexts
     * @param queryLookaheadBases look ahead this many bases during queries that produce cache misses
     */
    public MultiVariantDataSource(final List<FeatureInput<VariantContext>> featureInputs, final int queryLookaheadBases) {
        this(featureInputs, queryLookaheadBases, 0, 0, null, false);
    }

    /**
     * Creates a MultiVariantDataSource backed by the provided FeatureInputs. We will look ahead the specified number of bases
     * during queries that produce cache misses.
     *
     * @param featureInputs List of FeatureInput<VariantContext>> specifying sources of VariantContexts
     * @param queryLookaheadBases look ahead this many bases during queries that produce cache misses
     * @param cloudPrefetchBuffer  MB size of caching/prefetching wrapper for the data, if on Google Cloud (0 to disable).
     * @param cloudIndexPrefetchBuffer MB size of caching/prefetching wrapper for the index, if on Google Cloud (0 to disable).
     * @param reference reference to use when creating FeatureDataSources, may be null, only needed by GenomicsDB
     */
    public MultiVariantDataSource(final List<FeatureInput<VariantContext>> featureInputs, final int queryLookaheadBases, final int cloudPrefetchBuffer, final int cloudIndexPrefetchBuffer, final Path reference, final boolean skipDictionaryValidation) {
        Utils.validateArg(queryLookaheadBases >= 0, "Query lookahead bases must be >= 0");
        Utils.validateArg(featureInputs != null && featureInputs.size() > 0, "FeatureInputs list must be non-null and non-empty");

        featureInputs.forEach(
                featureInput -> featureDataSources.add(
                        new FeatureDataSource<>(featureInput, queryLookaheadBases, VariantContext.class, cloudPrefetchBuffer, cloudIndexPrefetchBuffer,
                                                reference)));

        // Ensure that the merged header and sequence dictionary that we use are in sync with each
        // other, and reflect the actual dictionaries used to do validation:
        //
        // 1) Cross validate the sequence dictionaries from each data source (which may be derived from an index
        //    in the case where its not embedded in the input file) to ensure they're mutually compatible
        // 2) Create and cache a merged header using versions of the individual headers from each data source that
        //    have been updated to include the actual dictionary returned from that data source
        //
        if (!skipDictionaryValidation) {
            validateAllSequenceDictionaries();
        }
        mergedHeader = getMergedHeader();
        mergedSamples = getSortedSamples();
        if ((mergedHeader == null || mergedHeader.getSequenceDictionary() == null) && featureInputs.size() > 1) {
            throw new UserException(
                    "No sequence dictionary was found for any input. When using multiple inputs, at least one input " +
                    "must have a sequence dictionary, or an index from which a sequence dictionary can be derived.");
        }

    }

    /**
     * Returns the aggregate sequence dictionary for this source of Variants. Uses the dictionary resulting
     * from merging available individual VCF headers (if present) for variant inputs.
     *
     * @return the sequence dictionary derived from the input sources, or null if no dictionary could be created
     * from any header or index file.
     */
    public SAMSequenceDictionary getSequenceDictionary() {
        return mergedHeader == null ? null : mergedHeader.getSequenceDictionary();
    }

    /**
     * Gets the merged header associated with this data source
     *
     * @return header associated with this data source as an Object
     */
    public VCFHeader getHeader() {
        return mergedHeader;
    }

    /**
     * Restricts traversals of this data source via {@link #iterator} to only return variants that overlap the provided
     * intervals. Calls to {@link #query(SimpleInterval)} are not affected by these intervals.
     *
     * Intervals MUST be non-overlapping and sorted in order of increasing start position, otherwise traversal
     * results will be incorrect.
     *
     * Passing in a null or empty interval List clears the intervals for traversal, making future iterations
     * over this data source unrestricted by intervals.
     *
     * @param intervals Our next full traversal will return only variants overlapping these intervals
     */
    public void setIntervalsForTraversal( final List<SimpleInterval> intervals ) {
        featureDataSources.forEach(ds -> ds.setIntervalsForTraversal(intervals));
    }

    /**
     * Gets an iterator over all variants in this data source, restricting traversal to variants
     * overlapping our intervals if intervals were provided via {@link #setIntervalsForTraversal(List)}
     *
     * Calling this method invalidates (closes) any previous iterator obtained from this method.
     *
     * @return an iterator over all variants in this data source, limited to variants that overlap the intervals supplied via {@link #setIntervalsForTraversal(List)} (if intervals were provided)
     */
    @Override
    public Iterator<VariantContext> iterator() {
        return getMergedIteratorFromDataSources(ds -> ds.iterator());
    }

    /**
     * Gets an iterator over all Variants in this data source that overlap the provided interval.
     *
     * This operation is not affected by intervals provided via {@link #setIntervalsForTraversal(List)}.
     *
     * Requires the backing files to have been indexed using the IndexFeatureFile tool, and to
     * be sorted in increasing order of start position for each contig.
     *
     * Query results are cached to improve the performance of future queries during typical access
     * patterns. See notes to the class as a whole for a description of the caching strategy.
     *
     * Calling this method potentially invalidates (closes) any other open iterator obtained
     * from this data source via a call to {@link #iterator}
     *
     * @param interval retrieve all Variants overlapping this interval
     * @return an iterator over all Variants in this data source that overlap the provided interval
     */
    @Override
    public Iterator<VariantContext> query( final SimpleInterval interval ) {
        return getMergedIteratorFromDataSources(ds -> ds.queryAndPrefetch(interval).iterator());
    }

    /**
     * Close any existing iterator, create a new iterator and update the local cached iterator reference.
     * @param iteratorFromSource function to retrieve individual iterator, to be applied to each data source
     * @return
     */
    private Iterator<VariantContext> getMergedIteratorFromDataSources(
            final Function<FeatureDataSource<VariantContext>, Iterator<VariantContext>> iteratorFromSource) {

        // Tribble documentation states that having multiple iterators open simultaneously over the same FeatureReader
        // results in undefined behavior
        closeOpenIterationIfNecessary();

        if (featureDataSources.size() > 1) {
            final List<CloseableIterator<VariantContext>> iterators = new ArrayList<>(featureDataSources.size());
            featureDataSources.forEach(ds -> iterators.add(getCloseableIteratorWrapper(iteratorFromSource.apply((ds)))));

            final VariantContextComparator varComparator = new VariantContextComparator(getSequenceDictionary());
            currentIterator = new MergingIterator<>(varComparator, iterators);
        } else {
            currentIterator = getCloseableIteratorWrapper(iteratorFromSource.apply(featureDataSources.get(0)));
        }
        return currentIterator;
    }

    /**
     * Get the logical name of this data source.
     *
     * @return the logical name of this data source
     */
    public String getName() {
        return "MultiVariantDataSource: ("
                + Utils.join(", ", featureDataSources.stream().map(fds -> fds.getName()).collect(Collectors.toList()))
                + ")";
    }

    /**
     * Permanently close this data source, invalidating any open iteration over it, and making it invalid for future
     * iterations and queries.
     */
    @Override
    public void close() {
        closeOpenIterationIfNecessary();
        featureDataSources.forEach(dataSource -> dataSource.close());
    }

    private SAMSequenceDictionary getMergedSequenceDictionary(VCFHeader header) {
        return header != null ? header.getSequenceDictionary() : null;
    }

    /**
     * Merge and sort the samples from each header requiring unique samples
     */
    private SortedSet<String> getSortedSamples() {
        final Map<String, VCFHeader> headers = featureDataSources
                .stream()
                .collect(Collectors.toMap(ds -> ds.getName(), ds -> (VCFHeader) ds.getHeader()));
        
        return VcfUtils.getSortedSampleSet(headers, GATKVariantContextUtils.GenotypeMergeType.REQUIRE_UNIQUE);
    }

    /**
     * Update each individual header with the sequence dictionary returned by the corresponding data source;
     * then merge the resulting headers.
     */
    private VCFHeader getMergedHeader() {
        final List<VCFHeader> headers = featureDataSources
                .stream()
                .map(ds -> getHeaderWithUpdatedSequenceDictionary(ds))
                .collect(Collectors.toList());

        // Now merge the headers using htsjdk, which is pretty promiscuous, and which only works properly
        // because of the cross-dictionary validation done in validateAllSequenceDictionaries.
        return headers.size() > 1 ?
                new VCFHeader(VCFUtils.smartMergeHeaders(headers, true)) :
                headers.get(0);
    }

    /**
     * We want the headers that are used to create the merged header to have the same sequence dictionary
     * that was  returned from the data source and used during validation (which may or may not be the
     * one that was embedded in the input file itself), so get the embedded one from the data source and
     * update it to include the actual sequence dictionary.
     */
    private VCFHeader getHeaderWithUpdatedSequenceDictionary(final FeatureDataSource<VariantContext> dataSource) {
        final VCFHeader header = (VCFHeader) dataSource.getHeader();
        final SAMSequenceDictionary sourceDict = dataSource.getSequenceDictionary();
        if (header.getSequenceDictionary() == null && sourceDict != null) {
            header.setSequenceDictionary(sourceDict);
        }
        return header;
    }

    /**
     * GATKTool only validates individual feature dictionaries against the reference dictionary, so cross-validate
     * all of the dictionaries against each other here by ensuring that each contig found in any dictionary has the
     * same length (and md5, when a value is present for that contig in both dictionaries) in every other dictionary
     * in which its present.
     */
    private void validateAllSequenceDictionaries() {
        final Map<String, FeatureDataSource<VariantContext>> contigMap = new HashMap<>();
        featureDataSources.forEach(
            ds -> getDataSourceDictionaryAndValidate(ds, contigMap)
        );
    }

    private void getDataSourceDictionaryAndValidate(final FeatureDataSource<VariantContext> ds,
                                                    final Map<String, FeatureDataSource<VariantContext>> contigMap) {
        final SAMSequenceDictionary dictionary = ds.getSequenceDictionary();
        if (dictionary == null) {
            logger.warn(
                    "A sequence dictionary is required for each input when using multiple inputs, and one could" +
                            " not be obtained for feature input: " + ds.getName() +
                            ". The input may not exist or may not have a valid header");
        } else {
            //This is HORRIFICALLY inefficient -- for tools with many inputs instead skip cross-validation by
            // overloading doDictionaryCrossValidation as false and require a reference
            dictionary.getSequences().forEach(
                    sourceSequence -> validateContigAgainstPreviousDataSource(sourceSequence, dictionary, contigMap, ds)
            );
        }
    }

    private void validateContigAgainstPreviousDataSource(final SAMSequenceRecord sourceSequence,
                  final SAMSequenceDictionary dictionary,
                  final Map<String, FeatureDataSource<VariantContext>> contigMap,
                  final FeatureDataSource<VariantContext> ds){
        final String sourceSequenceName = sourceSequence.getSequenceName();
        final FeatureDataSource<VariantContext> previousDataSource = contigMap.getOrDefault(sourceSequenceName, null);
        if (previousDataSource != null) {
            final SAMSequenceDictionary previousDictionary = previousDataSource.getSequenceDictionary();
            final SAMSequenceRecord previousSequence = previousDictionary.getSequence(sourceSequenceName);
            validateSequenceDictionaryRecords(
                    ds.getName(), dictionary, sourceSequence,
                    previousDataSource.getName(), previousDictionary, previousSequence);
        } else {
            contigMap.put(sourceSequenceName, ds);
        }
    }

    // Cross validate the length and md5 for a pair of sequence records.
    private void validateSequenceDictionaryRecords(
            final String sourceDataSourceName,
            final SAMSequenceDictionary sourceDictionary,
            final SAMSequenceRecord sourceSequence,
            final String targetDataSourceName,
            final SAMSequenceDictionary targetDictionary,
            final SAMSequenceRecord targetSequence)
    {
        if (!SequenceDictionaryUtils.sequenceRecordsAreEquivalent(sourceSequence, targetSequence)) {
            final String msg = String.format("Incompatible sequences found (%s: %d) and (%s: %d)",
                    sourceSequence.getSequenceName(),
                    sourceSequence.getSequenceLength(),
                    targetSequence.getSequenceName(),
                    targetSequence.getSequenceLength());
            throw new UserException.IncompatibleSequenceDictionaries(
                    msg,
                    sourceDataSourceName,
                    sourceDictionary,
                    targetDataSourceName,
                    targetDictionary
            );
        } else {
            final String targetMd5 = targetSequence.getMd5();
            final String sourceMd5 = sourceSequence.getMd5();
            if (Utils.xor(targetMd5 == null, sourceMd5 == null)) {
                final String msg = String.format(
                        "The MD5 value (%s) for sequence (%s) is present in at least one input sequence dictionary" +
                        " but missing in at least one other input sequence dictionary. In the case where an input VCF" +
                        " contains no embedded sequence dictionary, one may have been derived from the accompanying" +
                        " index; such derived dictionaries do not contain MD5 values and can result in this warning" +
                        " message.",
                        targetMd5 == null ? sourceMd5 : targetMd5,
                        sourceSequence.getSequenceName()
                );
                logger.warn(msg);
            } else if (!Objects.equals(targetMd5, sourceMd5)) {
                final String msg = String.format("Incompatible sequence MD5 values found (%s: %s) and (%s: %s)",
                        sourceSequence.getSequenceName(),
                        sourceMd5,
                        targetSequence.getSequenceName(),
                        targetMd5);
                throw new UserException.IncompatibleSequenceDictionaries(
                        msg,
                        sourceDataSourceName,
                        sourceDictionary,
                        targetDataSourceName,
                        targetDictionary
                );
            }
        }
    }

    /**
     * Close the iterator currently open over the data sources, if there is one.
     */
    private void closeOpenIterationIfNecessary() {
        if (currentIterator != null) {
            currentIterator.close();
            currentIterator = null;
        }
    }

    /**
     * Wrap the sourceIterator in a CloseableIterator to make it usable as a MergingIterator source.
     */
    private CloseableIterator<VariantContext> getCloseableIteratorWrapper(final Iterator<VariantContext> sourceIterator) {
        Utils.nonNull(sourceIterator);

        return new CloseableIterator<VariantContext>() {
            Iterator<VariantContext> delegateIterator = sourceIterator;
            @Override
            public void close() { delegateIterator = null; }

            @Override
            public boolean hasNext() {
                return delegateIterator != null && delegateIterator.hasNext();
            }

            @Override
            public VariantContext next() {
                if (!hasNext()) {
                    throw new NoSuchElementException("hasNext should be called before next");
                }
                return delegateIterator.next();
            }
        };
    }

    /**
     * Return lexicographically sorted set of uniquified sample names merged from across input data sources
     */
    public SortedSet<String> getSamples() {
        return mergedSamples;
    }
}
