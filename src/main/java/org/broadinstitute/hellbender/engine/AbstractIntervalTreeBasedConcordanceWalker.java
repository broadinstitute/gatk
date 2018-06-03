package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Locatable;
import htsjdk.tribble.Feature;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextComparator;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.commons.collections4.Predicate;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVIntervalTree;
import org.broadinstitute.hellbender.tools.walkers.validation.ConcordanceState;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Objects;

/**
 * A different implementation of variant call set concordance testing
 * (different from {@link AbstractConcordanceWalker}).
 *
 * This one is more appropriate for variants that have larger reference spans,
 * e.g. for CNV and SV.
 *
 * Essentially, the truth call set are traversed first to build an interval tree,
 * which is going to be tested for overlap later when traversing the evaluation call set.
 */
public abstract class AbstractIntervalTreeBasedConcordanceWalker extends GATKTool {

    public static final String TRUTH_VARIANTS_SHORT_NAME = "T";
    public static final String TRUTH_VARIANTS_LONG_NAME = "truth";

    public static final String EVAL_VARIANTS_SHORT_NAME = "E";
    public static final String EVAL_VARIANTS_LONG_NAME = "evaluation";

    public static final String CONFIDENCE_REGION_SHORT_NAME = "C";
    public static final String CONFIDENCE_REGION_LONG_NAME = "confidence";

    // The distance in bases to look ahead and cache when querying feature sources.
    public static final int CACHE_LOOKAHEAD = 100_000;

    public static final String FEATURE_NAME_FOR_TRUTH = "truth";

    // ARGUMENTS =======================================================================================================

    @Argument(doc = "TO BE IMPLEMENTED",
            shortName = CONFIDENCE_REGION_SHORT_NAME,
            fullName= CONFIDENCE_REGION_LONG_NAME,
            optional = true)
    protected String highConfidenceRegion;

    @Argument(doc = "path to a file holding mask on region to be excluded",
            shortName = "mask",
            fullName = "ref-region-mask-file",
            optional = true)
    protected String refRegionMaskFilePath;

    @Argument(doc = "A VCF containing truth variants",
            shortName = TRUTH_VARIANTS_LONG_NAME,
            fullName = TRUTH_VARIANTS_LONG_NAME)
    public String truthVariantsFile;

    @Argument(doc = "One VCF file containing call sets to be evaluated",
            shortName = EVAL_VARIANTS_SHORT_NAME,
            fullName = EVAL_VARIANTS_LONG_NAME)
    public String evalVariantsFile;

    // FIELDS ==========================================================================================================

    protected SVIntervalTree<String> maskedOut;
    protected SAMSequenceDictionary refSeqDict;
    private VariantContextComparator variantContextComparator;

    private FeatureDataSource<VariantContext> truthVariants;
    protected SVIntervalTree<List<TruthVariant>> truthTree;

    private FeatureDataSource<VariantContext> evalVariants;

    /**
     * For testing two objects of type T, using information stored in object of type R.
     * See example in {@link #getConcordanceTester()}
     */
    @FunctionalInterface
    public interface TriPredicate<T, R> {

        Boolean test(T t1, T t2, R r);

        default TriPredicate<T, R> and(final TriPredicate<? super T, ? super R> other) {
            Objects.requireNonNull(other);
            return (T t1, T t2, R r) -> test(t1, t2, r) && other.test(t1, t2, r);
        }

        default TriPredicate<T, R> negate() {
            return (T t1, T t2, R r) -> !test(t1, t2, r);
        }

        default TriPredicate<T, R> or(final TriPredicate<? super T, ? super R> other) {
            Objects.requireNonNull(other);
            return (T t1, T t2, R r) -> test(t1, t2, r) || other.test(t1, t2, r);
        }
    }

    protected TriPredicate<VariantContext, ReferenceContext> variantConcordanceTester;

    // TRIVIA ==========================================================================================================

    @Override
    public final SAMSequenceDictionary getBestAvailableSequenceDictionary() {
        return refSeqDict;
    }

    public final VCFHeader getTruthHeader() { return getHeader(truthVariants); }

    public final VCFHeader getEvalHeader() {
        return getHeader(evalVariants);
    }

    private static VCFHeader getHeader(final FeatureDataSource<VariantContext> source) {
        final Object header = source.getHeader();
        if ( ! (header instanceof VCFHeader) ) {
            throw new GATKException("Header for " + source.getName() + " is not in VCF header format");
        }
        return (VCFHeader)header;
    }

    /**
     * Convert from 1-based closed locatable [1, e]
     *         to 0-based semi-open interval [0, e).
     */
    protected final SVInterval convertLocatable(final Locatable locatable, final SAMSequenceDictionary sequenceDictionary) {
        return new SVInterval(sequenceDictionary.getSequenceIndex(locatable.getContig()),
                              locatable.getStart() - 1, locatable.getEnd());
    }

    // INITIALIZATION, TRAVERSAL, AND TERMINATION ======================================================================

    /**
     * A bunch of initializations.
     * The initializers can be overridden.
     */
    @Override
    protected final void onStartup() {

        IOUtils.assertFileIsReadable(Paths.get(truthVariantsFile));
        IOUtils.assertFileIsReadable(Paths.get(evalVariantsFile));

        super.onStartup();

        features = new FeatureManager(this, FeatureDataSource.DEFAULT_QUERY_LOOKAHEAD_BASES,
                cloudPrefetchBuffer, cloudIndexPrefetchBuffer, referenceArguments.getReferencePath());

        initializeMaskOrHighConfidence();

        initializeTruth();

        initializeEvalVariants();

        refSeqDict = truthVariants.getSequenceDictionary();
        variantContextComparator = getVariantContextComparator(refSeqDict);
        variantConcordanceTester = getConcordanceTester();
    }

    private void initializeMaskOrHighConfidence() {
        // sorry, I cannot find the BED file reader.....
        maskedOut = new SVIntervalTree<>();
        if (refRegionMaskFilePath != null) {
            try ( final BufferedReader rdr =
                          new BufferedReader(new InputStreamReader(BucketUtils.openFile(refRegionMaskFilePath))) ) {
                String line;
                while ( (line = rdr.readLine()) != null ) {
                    String[] bedLine = line.split("\t");
                    int sequenceIndex = refSeqDict.getSequenceIndex(bedLine[0]);
                    final int start = Integer.valueOf(bedLine[1]);
                    final int end = Integer.valueOf(bedLine[2]);
                    maskedOut.put(new SVInterval(sequenceIndex, start, end), "");
                }
            } catch (final IOException ioe ) {
                throw new GATKException("Unable to read intervals from " + refRegionMaskFilePath, ioe);
            }
        }
    }

    private void initializeTruth() {
        truthVariants = new FeatureDataSource<>(new FeatureInput<>(truthVariantsFile, FEATURE_NAME_FOR_TRUTH), CACHE_LOOKAHEAD, VariantContext.class);
        if (truthVariants.getSequenceDictionary() == null) {
            throw new UserException("Given truth variant VCF doesn't seem to have a valid reference sequence dictionary: " + truthVariantsFile);
        }
        if ( hasUserSuppliedIntervals() ) {
            truthVariants.setIntervalsForTraversal(userIntervals);
        }
        truthTree = buildIntervalTreeFromTruth();
    }

    private void initializeEvalVariants() {
        evalVariants = new FeatureDataSource<>(new FeatureInput<>(evalVariantsFile, FEATURE_NAME_FOR_TRUTH), CACHE_LOOKAHEAD, VariantContext.class);
        if ( hasUserSuppliedIntervals() ) {
            evalVariants.setIntervalsForTraversal(userIntervals);
        }

//        // when there are multiple eval call sets, one can do the following, and make changes accordingly
//        final List<FeatureInput<VariantContext>> evalVariantsFeatureInputs = new ArrayList<>(2);
//        final FeatureInput<VariantContext> featureInput = new FeatureInput<>(evalVariantsFile);
//        evalVariantsFeatureInputs.add(featureInput);
//        features.addToFeatureSources(0, featureInput, VariantContext.class, cloudPrefetchBuffer, cloudIndexPrefetchBuffer,
//                referenceArguments.getReferencePath());
//        MultiVariantDataSource multiVariantDataSource = new MultiVariantDataSource(evalVariantsFeatureInputs, CACHE_LOOKAHEAD, cloudPrefetchBuffer, cloudIndexPrefetchBuffer,
//                referenceArguments.getReferencePath());
    }

    /**
     * Evaluate each variant.
     *
     * One can override {@link #onTraversalStart()} for additional initialization,
     * and/or {@link #onTraversalSuccess()} for additional post-traversal tasks.
     */
    @Override
    public final void traverse() {
        int filteredEvalCount = 0;
        final Predicate<VariantContext> evalFilter = makeEvalVariantFilter();
        for (final VariantContext evalVariant : evalVariants) {
            if ( ! evalFilter.evaluate(evalVariant)) {++filteredEvalCount;continue;}

            final SimpleInterval variantInterval = new SimpleInterval(evalVariant);

            apply(evalVariant, null, getPaddedReferenceContext(variantInterval));
        }
        logger.info(filteredEvalCount + " eval variants skipped evaluation");
    }

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

    // OVERRIDE HIGHLY DESIRABLE =======================================================================================

    protected abstract void apply(final VariantContext eval, final ReadsContext readsContext, final ReferenceContext refContext);

    /**
     * Pad provided eval call accordingly for overlapping with the truth interval tree.
     */
    protected abstract SVInterval getPaddedSvInterval(final VariantContext eval);

    /**
     * Given {@code interval}, return padded reference context.
     */
    protected ReferenceContext getPaddedReferenceContext(final SimpleInterval interval) {
        return new ReferenceContext(reference, interval);
    }

    /**
     * Tools could override how padding are added,
     * if any padding is desired.
     */
    protected SVIntervalTree<List<TruthVariant>> buildIntervalTreeFromTruth() {
        final SVIntervalTree<List<TruthVariant>> truthTree = new SVIntervalTree<>();
        final Predicate<VariantContext> variantContextPredicate = makeTruthVariantFilter();
        final SAMSequenceDictionary sequenceDictionary = truthVariants.getSequenceDictionary();
        truthVariants.forEach(truth -> {
            if (!variantContextPredicate.evaluate(truth)) return;
            final SVInterval svInterval = new SVInterval(sequenceDictionary.getSequenceIndex(truth.getContig()),
                    truth.getStart() - 1, truth.getEnd());
            SVIntervalTree.Entry<List<TruthVariant>> listEntry = truthTree.find(svInterval);
            // pre-classify every truth call as FN, that is, without matching eval call, later could be modified accordingly by apply(...)
            final List<TruthVariant> value;
            if (listEntry == null) {
                value = Arrays.asList(new FalseNegative(truth));
            } else {
                value = new ArrayList<>(listEntry.getValue()); // listEntry is immutable
                value.add(new FalseNegative(truth));
            }
            truthTree.put(svInterval, value);
        });
        return truthTree;
    }

    protected Predicate<VariantContext> makeTruthVariantFilter() {
        return vc -> !vc.isFiltered();
    }

    /**
     * Override to filter out variants that should not be considered.
     * An example would be variants that are below the size threshold of 50 bp for SV call set evaluation.
     */
    protected Predicate<VariantContext> makeEvalVariantFilter() { return vc -> true; }

    /**
     * For comparing variants that don't overlap
     */
    protected VariantContextComparator getVariantContextComparator(final SAMSequenceDictionary dict) {
        return new VariantContextComparator(dict);
    }

    /**
     * Override this to customize the degree of agreement required to call a true positive.
     * Sometimes, for example, we may want just a single alt allele to agree
     * and sometime we may require all alt alleles to match.
     */
    protected abstract TriPredicate<VariantContext, ReferenceContext> getConcordanceTester();

    /**
     * Computes a score of concordance between the eval and truth.
     * This is particularly useful for variants that spans a large interval, i.e. CNV and SV calls.
     */
    protected abstract float computeConcordanceScore(final VariantContext eval, final VariantContext truth, final ReferenceContext context);

    // CORE CONCEPT OF THIS WALKER =====================================================================================

    public interface ClassifiedVariant extends Feature {
        ConcordanceState getConcordanceState();
    }
    public interface TruthVariant extends ClassifiedVariant {
        VariantContext getTruth();
    }
    public interface EvalVariant extends ClassifiedVariant {
        VariantContext getEval();
    }

    public static final class FalseNegative implements TruthVariant {
        private final VariantContext truth;

        public FalseNegative(final VariantContext truth) {
            this.truth = truth;
        }

        @Override
        public VariantContext getTruth() {
            return truth;
        }
        @Override
        public String getContig() {
            return truth.getContig();
        }
        @Override
        public int getStart() {
            return truth.getStart();
        }
        @Override
        public int getEnd() {
            return truth.getEnd();
        }

        @Override
        public ConcordanceState getConcordanceState() {
            return ConcordanceState.FALSE_NEGATIVE;
        }
    }

    public static abstract class EvalOnlyVariant implements EvalVariant {
        private final VariantContext eval;

        EvalOnlyVariant(final VariantContext eval) {
            this.eval = eval;
        }

        @Override
        public VariantContext getEval() {
            return eval;
        }
        @Override
        public String getContig() {
            return eval.getContig();
        }
        @Override
        public int getStart() {
            return eval.getStart();
        }
        @Override
        public int getEnd() {
            return eval.getEnd();
        }
    }
    public static final class FalsePositive extends EvalOnlyVariant {
        private ReasonForFalsePositive reasonForFalsePositive;

        public FalsePositive(final VariantContext eval, final ReasonForFalsePositive reasonForFalsePositive) {
            super(eval);
            this.reasonForFalsePositive = reasonForFalsePositive;
        }

        public ConcordanceState getConcordanceState() {
            return ConcordanceState.FALSE_POSITIVE;
        }

        public ReasonForFalsePositive getReasonForFalsePositive() {
            return reasonForFalsePositive;
        }

        public enum ReasonForFalsePositive {
            NO_OVERLAPPING_TRUTH("NO_OVP"),
            OVERLAPPING_TRUTH_WITH_DIFF_TYPE_OR_LOW_SUPPORT("WR_TYPE_OR_LOW_SUPP");

            private final String abbreviation;

            ReasonForFalsePositive(final String abbreviation) {
                this.abbreviation = abbreviation;
            }

            public String getAbbreviation() { return abbreviation; }
        }
    }
    public static final class FilteredTrueNegative extends EvalOnlyVariant {
        public FilteredTrueNegative(final VariantContext eval) {
            super(eval);
        }

        public ConcordanceState getConcordanceState() {
            return ConcordanceState.FILTERED_TRUE_NEGATIVE;
        }
    }

    public static abstract class BiVariant implements TruthVariant, EvalVariant {
        private final VariantContext eval;
        private final VariantContext supportingTruth;
        private final float maxSupportScore;

        BiVariant(final VariantContext eval, final VariantContext supportingTruth, final float maxSupportScore) {
            this.eval = eval;
            this.supportingTruth = supportingTruth;
            this.maxSupportScore = maxSupportScore;
        }

        @Override
        public final VariantContext getTruth() {
            return supportingTruth;
        }
        @Override
        public final VariantContext getEval() {
            return eval;
        }
        public final float getSupportScore() {
            return maxSupportScore;
        }

        @Override
        public final String getContig() {
            return eval.getContig();
        }
        @Override
        public final int getStart() {
            return eval.getStart();
        }
        @Override
        public final int getEnd() {
            return eval.getEnd();
        }
    }
    public static final class TruePositive extends BiVariant {

        public TruePositive(final VariantContext eval, final VariantContext supportingTruth, final float maxSupportScore) {
            super(eval, supportingTruth, maxSupportScore);
        }

        public final ConcordanceState getConcordanceState() {
            return ConcordanceState.TRUE_POSITIVE;
        }
    }
    public static final class FilteredFalseNegative extends BiVariant {

        public FilteredFalseNegative(final VariantContext eval, final VariantContext supportingTruth, final float maxSupportScore) {
            super(eval, supportingTruth, maxSupportScore);
        }

        public final ConcordanceState getConcordanceState() {
            return ConcordanceState.FILTERED_FALSE_NEGATIVE;
        }
    }
}
