package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.PeekableIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextComparator;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.commons.collections4.Predicate;
import org.apache.commons.collections4.iterators.FilterIterator;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;

import java.io.File;
import java.util.ArrayDeque;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Queue;
import java.util.Set;
import java.util.Spliterator;
import java.util.Spliterators;
import java.util.function.BiFunction;
import java.util.stream.StreamSupport;

/**
 * Base class for concordance walkers, which process one variant at a time from one or more sources of variants,
 * with optional contextual information from a reference, sets of reads, and/or supplementary sources of
 * Features.
 *
 * Subclasses must implement the {@link #apply} and {@link #makeMatchVariantsFunction} methods to process each variant
 * and may optionally implement {@link #onTraversalStart} and/or {@link #onTraversalSuccess}.
 *
 * Created by Takuto Sato 1/30/17, abstractified by David Benjamin on 2/22/17.
 * {@link #onTraversalStart}, {@link #onTraversalSuccess} and/or {@link #closeTool}.
 */
public abstract class AbstractConcordanceWalker extends WalkerBase {

    public static final String TRUTH_VARIANTS_LONG_NAME = "truth";
    public static final String EVAL_VARIANTS_SHORT_NAME = "eval";
    public static final String EVAL_VARIANTS_LONG_NAME = "evaluation";

    public static final String CONFIDENCE_REGION_LONG_NAME = "confidence";
    public static final String CONFIDENCE_REGION_SHORT_NAME = "C";
    
    // The distance in bases to look ahead and cache when querying feature sources.
    public static final int CACHE_LOOKAHEAD = 100_000;

    @Argument(shortName = TRUTH_VARIANTS_LONG_NAME, fullName = TRUTH_VARIANTS_LONG_NAME,
            doc = "A VCF containing truth variants", optional = false)
    public String truthVariantsFile;

    @Argument(shortName = EVAL_VARIANTS_SHORT_NAME, fullName = EVAL_VARIANTS_LONG_NAME,
            doc = "A VCF containing variants to be compared to the truth", optional = false)
    public String evalVariantsFile;

    @Argument(doc = "TO BE IMPLEMENTED",
            fullName= CONFIDENCE_REGION_LONG_NAME,
            shortName = CONFIDENCE_REGION_SHORT_NAME,
            optional = true)
    protected File highConfidenceRegion;

    private FeatureDataSource<VariantContext> truthVariants;
    private FeatureDataSource<VariantContext> evalVariants;
    protected SAMSequenceDictionary dict;
    private VariantContextComparator variantContextComparator;

    // Overriding the superclass method to favor the truth variants' dictionary
    @Override
    public final SAMSequenceDictionary getBestAvailableSequenceDictionary() {
        initializeTruthVariantsIfNecessary();
        return truthVariants.getSequenceDictionary();
    }

    private void initializeTruthVariantsIfNecessary() {
        if (! BucketUtils.fileExists(truthVariantsFile)) {
            throw new IllegalArgumentException("Truth variants file " + truthVariantsFile + " does not exist or is not readable.");
        }

        if (truthVariants == null) {
            truthVariants = new FeatureDataSource<>(new FeatureInput<>(truthVariantsFile, "truth"), CACHE_LOOKAHEAD, VariantContext.class, cloudPrefetchBuffer, cloudIndexPrefetchBuffer);
        }
    }


    protected Predicate<VariantContext> makeTruthVariantFilter() {
        return vc -> !vc.isFiltered();
    }

    // this may often be overridden
    protected Predicate<VariantContext> makeEvalVariantFilter() { return vc -> true; }

    /* this function takes a list of truth variant contexts (the first parameter) and a list of eval variant contexts (the second parameter) at the same positions, and returns a queue of TruthVersusEval objects which
    has matched the eval records with the truth records they should be compared to.  This comes into play if multiallelics have been split, for example.
     */
    protected BiFunction<List<VariantContext>, List<VariantContext>, Queue<TruthVersusEval>> makeMatchVariantsFunction() {
        return (u,v) -> {
            //match based onshouldVariantsBeMatched
            //we expect only a small number of variants in each list, so we loop through, as the overhead of a map removes its speed advantage
            final Queue<TruthVersusEval> queue = new ArrayDeque<>();
            Set<VariantContext> unmatchedTruth  = new HashSet<>(u);
            for (final VariantContext eval : v) {
                boolean matched = false;
                for (final VariantContext truth : u) {
                    if (shouldVariantsBeMatched(truth, eval)) {
                        if (!matched) {
                            queue.add(new TruthVersusEval(truth, eval));
                            matched = true;
                        }
                        unmatchedTruth.remove(truth);
                    }
                }
                if (!matched) {
                    queue.add(TruthVersusEval.evalOnly(eval));
                }
            }

            //deal with unmatched truth
            for (final VariantContext truth : unmatchedTruth) {
                queue.add(TruthVersusEval.truthOnly(truth));
            }

            return queue;
        };
    }

    protected abstract boolean shouldVariantsBeMatched(final VariantContext truth, final VariantContext eval);

    private Spliterator<TruthVersusEval> getSpliteratorForDrivingVariants() {
        final Iterator<VariantContext> truthIterator = new FilterIterator<>(truthVariants.iterator(), makeTruthVariantFilter());
        final Iterator<VariantContext> evalIterator = new FilterIterator<>(evalVariants.iterator(), makeEvalVariantFilter());
        return new ConcordanceIterator(truthIterator, evalIterator, makeMatchVariantsFunction()).spliterator();
    }

    // ********** The basic traversal structure of GATKTool
    @Override
    protected final void onStartup() {
        super.onStartup();

        initializeTruthVariantsIfNecessary();
        evalVariants = new FeatureDataSource<>(new FeatureInput<>(evalVariantsFile, "eval"), CACHE_LOOKAHEAD, VariantContext.class, cloudPrefetchBuffer, cloudIndexPrefetchBuffer);

        if ( hasUserSuppliedIntervals() ) {
            truthVariants.setIntervalsForTraversal(userIntervals);
            evalVariants.setIntervalsForTraversal(userIntervals);
        }
        dict = getBestAvailableSequenceDictionary();
        variantContextComparator = new VariantContextComparator(dict);
    }

    // one can override onTraversalStart() for additional initialization

    // the primary work of the walker.  Must be overridden in implementing classes.
    protected abstract void apply(final TruthVersusEval truthVersusEval, final ReadsContext readsContext, final ReferenceContext refContext);

    /**
     * {@inheritDoc}
     *
     * Implementation of concordance traversal.
     *
     * NOTE: You should only override {@link #traverse()} if you are writing a new walker base class in the
     * engine package that extends this class. It is not meant to be overridden by tools outside of the
     * engine package.
     */
    @Override
    public void traverse() {
        // Process each variant in the input stream.
        StreamSupport.stream(getSpliteratorForDrivingVariants(), false)
                .forEach(truthVersusEval -> {
                    final SimpleInterval variantInterval = new SimpleInterval(truthVersusEval);
                    apply(truthVersusEval, new ReadsContext(reads, variantInterval), new ReferenceContext(reference, variantInterval));
                    progressMeter.update(variantInterval);
                });
    }

    // one can override onTraversalSuccess() for additional post-traversal tasks

    @Override
    protected final void onShutdown() {
        super.onShutdown();
        if( truthVariants != null ) {
            truthVariants.close();
        }
        if( evalVariants != null) {
            evalVariants.close();
        }
    }
    // ********** End of basic traversal methods

    public final VCFHeader getTruthHeader() { return getHeader(truthVariants); }

    public final VCFHeader getEvalHeader() { return getHeader(evalVariants); }

    private static VCFHeader getHeader(final FeatureDataSource<VariantContext> source) {
        final Object header = source.getHeader();
        if ( ! (header instanceof VCFHeader) ) {
            throw new GATKException("Header for " + source.getName() + " is not in VCF header format");
        }
        return (VCFHeader)header;
    }

    /**
     * encapsulate the iteration over truth and eval.  hasNext() returns true if either the truth or eval iterator have
     * any variants remaining.  next() produces variants in genomic order, along with their concordance state of
     * false positive for an eval variant without concordant truth variant, false negative for a truth variant without
     * a concordant eval variant, and true positive for concordant eval and truth variants.
     *
     * This implementation correctly consumes a single next() of the wrapped iterators for false positives and
     * false negatives and consumes a next() for both iterators in case of a true positive.
     */
    private class ConcordanceIterator implements Iterator<TruthVersusEval> {
        private final PeekableIterator<VariantContext> truthIterator;
        private final PeekableIterator<VariantContext> evalIterator;
        private final BiFunction<List<VariantContext>, List<VariantContext>, Queue<TruthVersusEval>> matchVariantsFunction;
        private final Queue<TruthVersusEval> truthVersusEvalQueue = new ArrayDeque<>();

        protected ConcordanceIterator(final  Iterator<VariantContext> truthIterator, final Iterator<VariantContext> evalIterator,
                                      final BiFunction<List<VariantContext>, List<VariantContext>, Queue<TruthVersusEval>> matchVariantsFunction) {
            this.truthIterator = new PeekableIterator<>(truthIterator);
            this.evalIterator = new PeekableIterator<>(evalIterator);
            this.matchVariantsFunction = matchVariantsFunction;
        }

        public boolean hasNext() { return truthIterator.hasNext() || evalIterator.hasNext() || !truthVersusEvalQueue.isEmpty(); }

        public TruthVersusEval next() {

            if (!truthVersusEvalQueue.isEmpty()) {
                return truthVersusEvalQueue.poll();
            }

            if (!truthIterator.hasNext()) {
                return TruthVersusEval.evalOnly(evalIterator.next());
            } else if (!evalIterator.hasNext()) {
                return TruthVersusEval.truthOnly(truthIterator.next());
            }

            final int positionCompare = variantContextComparator.compare(truthIterator.peek(), evalIterator.peek());
            if (positionCompare > 0) {
                return TruthVersusEval.evalOnly(evalIterator.next());
            } else if (positionCompare < 0) {
                return TruthVersusEval.truthOnly(truthIterator.next());
            } else {
                //get all eval records at this location
                final List<VariantContext> evalVariants = getNextVariantGroup(evalIterator);

                //get all truth records at this location
                final List<VariantContext> truthVariants = getNextVariantGroup(truthIterator);

                //if there is only one truth and eval at this location, can just return those, otherwise need to match using matching function
                if (evalVariants.size() == 1 && truthVariants.size() == 1) {
                    return new TruthVersusEval(evalVariants.get(0), truthVariants.get(0));
                } else {
                    truthVersusEvalQueue.addAll(matchVariantsFunction.apply(evalVariants, truthVariants));
                    return truthVersusEvalQueue.poll();
                }
            }
        }

        private List<VariantContext> getNextVariantGroup(final PeekableIterator<VariantContext> iterator) {
            final List<VariantContext> variants = new LinkedList<>();
            VariantContext currentVariant = iterator.next();
            variants.add(currentVariant);
            while(iterator.hasNext() && variantContextComparator.compare(currentVariant, iterator.peek()) == 0) {
                currentVariant = iterator.next();
                variants.add(currentVariant);
            }

            return variants;
        }

        private Spliterator<TruthVersusEval> spliterator() {
            return Spliterators.spliteratorUnknownSize(this, 0);
        }
    }

    /**
     * store a truth vc in case of a false negative, an eval vc in case of a false positive, or a concordance pair of
     * truth and eval in case of a true positive.
     */
    protected static class TruthVersusEval implements Locatable {
        protected final VariantContext truth;
        protected final VariantContext eval;

        public TruthVersusEval(final VariantContext truth, final VariantContext eval) {
            this.truth = truth;
            this.eval = eval;
        }

        public static TruthVersusEval evalOnly(final VariantContext eval) {
            return new TruthVersusEval(null, eval);
        }

        public static TruthVersusEval truthOnly(final VariantContext truth) {
            return new TruthVersusEval(truth, null);
        }

        public VariantContext getTruth() {
            Utils.validateArg(truth != null, () -> "This is variant has no truth VariantContext.");
            return truth;
        }

        public VariantContext getEval() {
            Utils.validateArg(eval != null, () -> "This is variant has no eval VariantContext.");
            return eval;
        }

        public boolean hasTruth() {return truth != null; }
        public boolean hasEval() { return eval != null; }

        public VariantContext getTruthIfPresentElseEval() { return truth != null ? truth : eval; }

        public String getContig() { return getTruthIfPresentElseEval().getContig(); }
        public int getStart() { return getTruthIfPresentElseEval().getStart(); }
        public int getEnd() { return getTruthIfPresentElseEval().getEnd(); }
    }

}
