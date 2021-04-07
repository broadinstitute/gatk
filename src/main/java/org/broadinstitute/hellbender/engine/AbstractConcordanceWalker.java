package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.PeekableIterator;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextComparator;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.commons.collections4.Predicate;
import org.apache.commons.collections4.iterators.FilterIterator;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.walkers.validation.ConcordanceState;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;

import java.io.File;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Optional;
import java.util.Spliterator;
import java.util.Spliterators;
import java.util.stream.StreamSupport;

/**
 * Base class for concordance walkers, which process one variant at a time from one or more sources of variants,
 * with optional contextual information from a reference, sets of reads, and/or supplementary sources of
 * Features.
 *
 * Subclasses must implement the {@link #apply} and {@link #areVariantsAtSameLocusConcordant} methods to process each variant
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

    private Spliterator<TruthVersusEval> getSpliteratorForDrivingVariants() {
        final Iterator<VariantContext> truthIterator = new FilterIterator<>(truthVariants.iterator(), makeTruthVariantFilter());
        final Iterator<VariantContext> evalIterator = new FilterIterator<>(evalVariants.iterator(), makeEvalVariantFilter());
        return new ConcordanceIterator(truthIterator, evalIterator).spliterator();
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
        final LinkedHashMap<List<Allele>, VariantContext> nextEvalVariants = new LinkedHashMap<>();

        protected ConcordanceIterator(final  Iterator<VariantContext> truthIterator, final Iterator<VariantContext> evalIterator) {
            this.truthIterator = new PeekableIterator<>(truthIterator);
            this.evalIterator = new PeekableIterator<>(evalIterator);
        }

        public boolean hasNext() { return truthIterator.hasNext() || evalIterator.hasNext() || !nextEvalVariants.isEmpty(); }

        public TruthVersusEval next() {

            if (!truthIterator.hasNext()) {
                if (!nextEvalVariants.isEmpty()) {
                    final VariantContext evalVariant = nextEvalVariants.values().iterator().next();
                    nextEvalVariants.remove(evalVariant.getAlleles());
                    return evalVariant.isFiltered()? TruthVersusEval.filteredTrueNegative(evalVariant) : TruthVersusEval.falsePositive(evalVariant);
                } else {
                    return nextEvalOnlyVariant();
                }
            } else if (!evalIterator.hasNext() && nextEvalVariants.isEmpty()) {
                return nextTruthOnlyVariant();
            }

            final int positionCompare = variantContextComparator.compare(truthIterator.peek(), nextEvalVariants.isEmpty()? evalIterator.peek() : nextEvalVariants.values().iterator().next());
            if (positionCompare > 0) {
                if (!nextEvalVariants.isEmpty()) {
                    final VariantContext evalVariant = nextEvalVariants.values().iterator().next();
                    nextEvalVariants.remove(evalVariant.getAlleles());
                    return evalVariant.isFiltered()? TruthVersusEval.filteredTrueNegative(evalVariant) : TruthVersusEval.falsePositive(evalVariant);
                } else {
                    return nextEvalOnlyVariant();
                }
            } else if (positionCompare < 0) {
                return nextTruthOnlyVariant();
            } else {
                if (nextEvalVariants.isEmpty()) {
                    VariantContext currentEvalVariant = evalIterator.next();
                    nextEvalVariants.put(currentEvalVariant.getAlleles(), currentEvalVariant);
                    while (evalIterator.hasNext() && variantContextComparator.compare(currentEvalVariant, evalIterator.peek()) == 0) {
                        currentEvalVariant = evalIterator.next();
                        nextEvalVariants.put(currentEvalVariant.getAlleles(), currentEvalVariant);
                    }
                }


                final VariantContext truthVariant = truthIterator.next();
                final VariantContext evalVariant = nextEvalVariants.remove(truthVariant.getAlleles());
                if (evalVariant == null) {
                    return TruthVersusEval.falseNegative(truthVariant);
                } else if (evalVariant.isFiltered()) {
                    return TruthVersusEval.filteredFalseNegative(truthVariant, evalVariant);
                } else if (areVariantsAtSameLocusConcordant(truthVariant, evalVariant)) {
                    return TruthVersusEval.truePositive(truthVariant, evalVariant);
                } else {
                    nextEvalVariants.put(evalVariant.getAlleles(), evalVariant);
                    return TruthVersusEval.falseNegative(truthVariant);
                }
            }
        }

        // if a variant only exists in the eval vcf, it is either a false positive or a filtered true negative
        private TruthVersusEval nextEvalOnlyVariant() {
            final VariantContext evalVariant = evalIterator.next();
            return evalVariant.isFiltered() ? TruthVersusEval.filteredTrueNegative(evalVariant) : TruthVersusEval.falsePositive(evalVariant);
        }

        private TruthVersusEval nextTruthOnlyVariant() { return TruthVersusEval.falseNegative(truthIterator.next()); }

        private Spliterator<TruthVersusEval> spliterator() {
            return Spliterators.spliteratorUnknownSize(this, 0);
        }
    }

    // override this to customize the degree of agreement required to call a true positive.  Sometimes, for example, we may want
    // just a single alt allele to agree and sometime we may require all alt alleles to match
    protected abstract boolean areVariantsAtSameLocusConcordant(final VariantContext truth, final VariantContext eval);

    /**
     * store a truth vc in case of a false negative, an eval vc in case of a false positive, or a concordance pair of
     * truth and eval in case of a true positive.
     */
    protected static class TruthVersusEval implements Locatable {
        private final Optional<VariantContext> truth;
        private final Optional<VariantContext> eval;
        private final ConcordanceState concordanceState;

        private TruthVersusEval(final Optional<VariantContext> truth, final Optional<VariantContext> eval,
                                final ConcordanceState concordanceState) {
            this.truth = truth;
            this.eval = eval;
            this.concordanceState = concordanceState;
        }

        public static TruthVersusEval falseNegative(final VariantContext truth) {
            return new TruthVersusEval(Optional.of(truth), Optional.empty(), ConcordanceState.FALSE_NEGATIVE);
        }

        public static TruthVersusEval falsePositive(final VariantContext eval) {
            return new TruthVersusEval(Optional.empty(), Optional.of(eval), ConcordanceState.FALSE_POSITIVE);
        }

        public static TruthVersusEval truePositive(final VariantContext truth, final VariantContext eval) {
            return new TruthVersusEval(Optional.of(truth), Optional.of(eval), ConcordanceState.TRUE_POSITIVE);
        }

        public static TruthVersusEval filteredFalseNegative(final VariantContext truth, final VariantContext eval) {
            return new TruthVersusEval(Optional.of(truth), Optional.of(eval), ConcordanceState.FILTERED_FALSE_NEGATIVE);
        }

        public static TruthVersusEval filteredTrueNegative(final VariantContext eval) {
            return new TruthVersusEval(Optional.empty(), Optional.of(eval), ConcordanceState.FILTERED_TRUE_NEGATIVE);
        }

        public ConcordanceState getConcordance() {
            return concordanceState;
        }

        public VariantContext getTruth() {
            Utils.validateArg(truth.isPresent(), () -> "This is a " + concordanceState.toString() + " and has no truth VariantContext.");
            return truth.get();
        }

        public VariantContext getEval() {
            Utils.validateArg(eval.isPresent(), () -> "This is a " + concordanceState.toString() + " and has no eval VariantContext.");
            return eval.get();
        }

        public boolean hasTruth() {return truth.isPresent(); }
        public boolean hasEval() { return eval.isPresent(); }

        public VariantContext getTruthIfPresentElseEval() { return truth.orElseGet(() -> eval.get()); }

        public String getContig() { return getTruthIfPresentElseEval().getContig(); }
        public int getStart() { return getTruthIfPresentElseEval().getStart(); }
        public int getEnd() { return getTruthIfPresentElseEval().getEnd(); }
    }

}
