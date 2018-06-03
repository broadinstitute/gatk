package org.broadinstitute.hellbender.tools.walkers.sv;

import com.google.common.collect.Sets;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.collections4.Predicate;
import org.apache.commons.lang.StringUtils;
import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.AbstractIntervalTreeBasedConcordanceWalker;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SimpleSVType;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVIntervalTree;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVUtils;
import org.broadinstitute.hellbender.tools.walkers.validation.ConcordanceState;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;

import java.io.BufferedOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.*;
import java.util.stream.Collectors;

/**
 * (Internal) Evaluates concordance between given call set and a single unique truth set for non-genotyped SV vcf
 *
 * <p>
 *     This is an experimental tool and should not be of interest to most researchers.
 *     It is a prototype of evaluating non-genotyped structural variant calls against a single truth,
 *     and has been under active developments.
 * </p>
 *
 * <p>
 *     This tool pads the interval around the interval bounded by POS and END as given in each record,
 *     by an amount specified in the argument "padding" to look for matching truth record, then
 *     test if the haplotypes by the truth and the evaluation calls are concordant.
 * </p>
 *
 * <p>
 *     The VCF files can be genotyped, but the FORMAT columns are currently ignored.
 * </p>
 *
 * <h3>Inputs</h3>
 * <ul>
 *     <li>path to single VCF file containing the truth</li>
 *     <li>path to single VCF file to be evaluated</li>
 *     <li>(optional) path to interval file containing high confidence region, and/or</li>
 *     <li>(optional) path to interval file containing regions to be excluded</li>
 *     <li>(optional) path to file containing non-canonical chromosomes (e.g. chromosomes other than chr1-22, chrX and chrY in hg38 reference)</li>
 *     <li></li>
 * </ul>
 *
 * <h3>Output (currently)</h3>
 * <ul>
 *     <li>prints to screen a summary table containing summary on the TP/FP/FN/FFN/FFP numbers.</li>
 *     <li>writes a BED-4 file containing: TP, FP, FILTERED_TRUE_NEGATIVE (record correctly filtered), and FILTERED_FALSE_NEGATIVE</li>
 * </ul>
 *
 * <h3>Usage example</h3>
 * <pre>
 *   gatk EvaluateConcordanceOnNonGenotypedSVVCFs \
 *     -T truth.vcf \
 *     -E evalCallSet.vcf \
 *     -O outputPrefix
 * </pre>
 */
@DocumentedFeature
@BetaFeature
@CommandLineProgramProperties(
        oneLineSummary = "(Internal) Evaluates concordance between given call set and a single unique truth set for non-genotyped SV vcf",
        summary = "This tool is used in development and should not be of interest to most researchers. " +
                "It is a prototype of evaluating non-genotyped structural variant calls against a single truth, " +
                "and has been under active developments. " +
                "This tool pads the interval around the interval bounded by POS and END as given in each record, " +
                "by an amount specified in the argument \"padding\" to look for matching truth record, then " +
                "tests if the haplotypes by the truth and the evaluation calls are concordant. ",
        programGroup = StructuralVariantDiscoveryProgramGroup.class)
public final class EvaluateConcordanceOnNonGenotypedSVVCFs extends AbstractIntervalTreeBasedConcordanceWalker {

    @Argument(doc = "file containing non-canonical chromosome names (e.g chrUn_KI270588v1) in the reference; " +
            "all chromosome will be considered canonical if not provided",
            shortName = "alt-tigs",
            fullName = "non-canonical-contig-names-file",
            optional = true)
    private String nonCanonicalChromosomeNamesFile;

    @Advanced
    @Argument(doc = "indicate evaluation is being performed on non-human calls",
            shortName = "nonhs",
            fullName = "non-human",
            optional = true)
    private Boolean nonHumanVCF = false;

    @Advanced
    @Argument(doc = "string value for chromosome Y; takes value of \"Y\" or \"chrY\", or when left unprovided, " +
            "evaluation is skipped on chromosome Y; applicable only when \"non-human\" is false",
            shortName = "chry",
            fullName = "chr-y-string",
            optional = true)
    private String chrYString;

    /**
     * A fixed amount of padding for each VCF record from the eval call set.
     * Each individual record might be padded (more or less bases) as judged appropriate by the tool itself.
     */
    @Argument(doc = "Padding around variants for overlapping",
            fullName = "padding",
            optional = true)
    private Integer fixedPadding = 50;

    @Argument(doc = "prefix for output",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME)
    private String outPrefix;

    @Argument(doc = "Minimum reciprocal overlap fraction for two concordant variants of the type DEL and INV",
            fullName = "recip-ovp-frac",
            optional = true)
    private Float recipOvpFrac = 0.5f;

    @Argument(doc = "Maximum difference in length between two concordant variants or type INS and/or DUP, " +
            "i.e. the difference in their length cannot exceed this fraction of the shorter of the two",
            fullName = "max-ins-len-diff-frac",
            optional = true)
    private Float insLenDiffFrac = 0.1f;

    @Advanced
    @Argument(doc = "Minimum mapping quality of the evidence contig mappings for a eval VCF record to be considered;" +
            " setting to zero turns off this filter",
            fullName = "min-map-qual",
            optional = true)
    private Integer minMQ = 0;

    @Advanced
    @Argument(doc = "Minimum alignment length of the evidence contig mappings; setting to zero turns off this filter",
            fullName = "min-aln-len",
            optional = true)
    private Integer minAlnLen = 0;

    // =================================================================================================================

    private static final Set<String> CHROMOSOME_Y_NAMES = Sets.newHashSet("Y", "chrY");
    private Set<String> chromosomesOfInterest = new HashSet<>();

    // TODO: 6/4/18 this is initialized in a way that $ 1 - insLenDiffScoreCoeff * lenDiff/Math.min(evalLen, truthLen) < recipOvpFrac $ for discordant insertion/duplication
    private Float insLenDiffScoreCoeff;

    private final List<EvalVariant> classifiedVariants = new ArrayList<>(50_000); // guess

    // =================================================================================================================

    @Override
    public final void onTraversalStart() {
        insLenDiffScoreCoeff = (1f - recipOvpFrac)/insLenDiffFrac;
        if (nonCanonicalChromosomeNamesFile != null) {
            chromosomesOfInterest = SVUtils.getCanonicalChromosomes(nonCanonicalChromosomeNamesFile, refSeqDict);
        } else {
            chromosomesOfInterest = refSeqDict.getSequences().stream().map(SAMSequenceRecord::getSequenceName)
                                              .collect(Collectors.toCollection(HashSet::new));
        }
        if (!nonHumanVCF && StringUtils.isEmpty(chrYString)) {
            chromosomesOfInterest.removeAll(CHROMOSOME_Y_NAMES);
        }
    }

    @Override
    protected Predicate<VariantContext> makeEvalVariantFilter() {
        return vc -> {
            // currently only considers variants on canonical chromosomes, filter turned off if chromosomesOfInterest is all chromosomes in VCF header
            if ( ! chromosomesOfInterest.contains(vc.getContig()) )
                return false;

            if ( maskedOut.hasOverlapper(convertLocatable(vc, refSeqDict)) )  // masked out
                return false;

            if ( vc.hasAttribute(GATKSVVCFConstants.MAPPING_QUALITIES) ) {
                int maxMQ = SVUtils.getAttributeAsStringStream(vc, GATKSVVCFConstants.MAPPING_QUALITIES).mapToInt(Integer::new).max().orElse(0);
                if  (maxMQ < minMQ)
                    return false;
            }
            if ( vc.hasAttribute(GATKSVVCFConstants.MAX_ALIGN_LENGTH) ) {
                if (vc.getAttributeAsInt(GATKSVVCFConstants.MAX_ALIGN_LENGTH, 0) < minAlnLen)
                    return false;
            }

            final int svLen = Math.abs(vc.getAttributeAsInt(GATKSVVCFConstants.SVLEN, 0));
            final String svType = vc.getAttributeAsString(GATKSVVCFConstants.SVTYPE, "");
            if (svType.equals(SimpleSVType.SupportedType.DEL.name())) { // some tiny deletion calls are actually RPL calls with micro deletions
                if ( (!vc.hasAttribute(GATKSVVCFConstants.INSERTED_SEQUENCE)) && svLen < StructuralVariationDiscoveryArgumentCollection.STRUCTURAL_VARIANT_SIZE_LOWER_BOUND) // no ins_seq AND tiny
                    return false;
            } else if (svType.equals(SimpleSVType.SupportedType.INS.name()) || svType.equals(SimpleSVType.SupportedType.DUP.name())) {
                if (svLen < StructuralVariationDiscoveryArgumentCollection.STRUCTURAL_VARIANT_SIZE_LOWER_BOUND)
                    return false;
            }

            return true;
        };
    }

    /**
     * A 2-step evaluation mechanism:
     * <ul>
     *     <li> first test if {@code eval} (padded by {@link #getPaddedSvInterval(VariantContext)}) overlaps with any record in the truth tree </li>
     *     <li> further test concordance if there's such overlap </li>
     * </ul>
     */
    @Override
    protected void apply(final VariantContext eval, final ReadsContext readsContext, final ReferenceContext refContext) {
        final SVInterval evalInterval = getPaddedSvInterval(eval);
        if ( truthTree.hasOverlapper(evalInterval) ) {
            final SupportingTruth supportingTruth = getConcordantTruthCalls(eval, evalInterval);
            if (supportingTruth == null) {
                classifiedVariants.add(new FalsePositive(eval, FalsePositive.ReasonForFalsePositive.OVERLAPPING_TRUTH_WITH_DIFF_TYPE_OR_LOW_SUPPORT));
            } else {
                classifiedVariants.add( eval.isFiltered() ? new FilteredFalseNegative(eval, supportingTruth.truthVar, supportingTruth.maxSupportScore)
                                                          : new TruePositive(eval, supportingTruth.truthVar, supportingTruth.maxSupportScore));
            }
        } else {
            if (eval.isFiltered()) {
                classifiedVariants.add(new FilteredTrueNegative(eval));
            } else {
                classifiedVariants.add(new FalsePositive(eval, FalsePositive.ReasonForFalsePositive.NO_OVERLAPPING_TRUTH));
            }
        }
    }

    // TODO: 6/3/18 padding can be improved
    @Override
    protected SVInterval getPaddedSvInterval(final VariantContext eval) {

        final int homologyLength = eval.getAttributeAsInt(GATKSVVCFConstants.HOMOLOGY_LENGTH, 0);

        int insTypeRightPadding = 0;
        final String svType = eval.getAttributeAsString(GATKSVVCFConstants.SVTYPE, "");
        if (svType.equals(GATKSVVCFConstants.SYMB_ALT_ALLELE_INS) || svType.equals(GATKSVVCFConstants.SYMB_ALT_ALLELE_DUP))
            insTypeRightPadding = eval.getAttributeAsInt(GATKSVVCFConstants.SVLEN, 0);

        // we only pad to the right with homology because we expect variants to follow left-align convention (GATK-SV does this)
        final int leftPadding = fixedPadding;
        final int rightPadding = fixedPadding + homologyLength + insTypeRightPadding;
        final int chromosomeLength = refSeqDict.getSequence(eval.getContig()).getSequenceLength();
        return new SVInterval(refSeqDict.getSequenceIndex(eval.getContig()),
                              Math.max(0, eval.getStart() - leftPadding - 1),
                              Math.min(eval.getEnd() + rightPadding, chromosomeLength));
    }

    private static final class SupportingTruth {

        private final VariantContext truthVar;
        private final float maxSupportScore;

        SupportingTruth(final VariantContext truthVar, final float maxSupportScore) {
            this.truthVar = truthVar;
            this.maxSupportScore = maxSupportScore;
        }
    }

    /**
     * Given an eval variant and its corresponding padded interval, return concordance truth variant, if any.
     * @param eval          variant to be evaluated
     * @param evalInterval  padded interval constructed from {@code eval}
     * @return              {@code null} if no concordant truth call exist
     */
    private SupportingTruth getConcordantTruthCalls(final VariantContext eval, final SVInterval evalInterval) {
        SupportingTruth bestSupportingTruth = null;

        final Iterator<SVIntervalTree.Entry<List<TruthVariant>>> overlappingTruths = truthTree.overlappers(evalInterval);

        Float maxConcordanceScore = 0f;
        while (overlappingTruths.hasNext()) {
            final List<TruthVariant> truthVariants = overlappingTruths.next().getValue();
            for (int i = 0; i < truthVariants.size(); ++i ) {
                final TruthVariant truthVariant = truthVariants.get(i);
                final VariantContext truthVarCtx = truthVariant.getTruth();
                float concordanceScore = computeConcordanceScore(eval, truthVarCtx, null);
                final Boolean concordant = variantConcordanceTester.test(eval, truthVarCtx, null);
                if (concordant && concordanceScore > maxConcordanceScore) {
                    SupportingTruth supportingTruth = new SupportingTruth(truthVarCtx, concordanceScore);
                    final TruthVariant updatedTruthVariant =  eval.isFiltered() ? new FilteredFalseNegative(eval, supportingTruth.truthVar, supportingTruth.maxSupportScore)
                                                                                : new TruePositive(eval, supportingTruth.truthVar, supportingTruth.maxSupportScore);
                    truthVariants.set(i, updatedTruthVariant);
                    bestSupportingTruth = supportingTruth;
                }
            }
        }
        return bestSupportingTruth;
    }

    // TODO: 6/3/18 test if heuristic overlap-based concordance score is over a hard threshold, in the long run, haplotype concordance should be evaluated instead
    /**
     * Assuming argument VariantContext's overlap.
     */
    @Override
    protected TriPredicate<VariantContext, ReferenceContext> getConcordanceTester() {
        return (eval, truth, r) -> recipOvpFrac <= computeConcordanceScore(eval, truth, r);
    }

    /**
     * >= 0.5 is considered concordant.
     */
    @Override
    protected float computeConcordanceScore(final VariantContext eval, final VariantContext truth, final ReferenceContext context) {
        final String evalType = eval.getAttributeAsString(GATKSVVCFConstants.SVTYPE, "");
        final String truthType = truth.getAttributeAsString(GATKSVVCFConstants.SVTYPE, "");

        // TODO: 6/3/18 we only consider insertion length but DO NOT YET MATCH inserted sequence (we should, in the long run, but that involves haplotype reconstruction)
        if (evalType.equals(SimpleSVType.SupportedType.INS.name()) || evalType.equals(SimpleSVType.SupportedType.DUP.name())) {
            if (truthType.equals(SimpleSVType.SupportedType.INS.name())) { // TODO: 6/3/18 truth type doesn't call DUP, instead type them as INS
                int evalLen = eval.getAttributeAsInt(GATKSVVCFConstants.SVLEN, 0);
                int truthLen = truth.getAttributeAsInt(GATKSVVCFConstants.SVLEN, 0);
                float diff = Math.abs(evalLen - truthLen);
                return 1 - insLenDiffScoreCoeff * Math.max(diff/evalLen, diff/truthLen);
            } else
                return 0f;
        } else if (evalType.equals(SimpleSVType.SupportedType.DEL.name())){

            if (truthType.equals(SimpleSVType.SupportedType.DEL.name())) { // NOTE: one can use SVLEN, but remember the truth set might report positive SVLEN for DEL
                SVInterval evalRange = new SVInterval(0, eval.getStart() - 1, eval.getEnd());
                SVInterval truthRange = new SVInterval(0, truth.getStart() - 1, truth.getEnd());
                float overlapLen = evalRange.overlapLen(truthRange);
                return Math.min(overlapLen/evalRange.getLength(), overlapLen/truthRange.getLength());
            } else if (truthType.equals(SimpleSVType.SupportedType.INS.name()) && eval.hasAttribute(GATKSVVCFConstants.INSERTED_SEQUENCE_LENGTH)) { // some call set might report replacement as deletion with insertion
                final int truthInsLen = truth.getAttributeAsInt(GATKSVVCFConstants.SVLEN, 0);
                final int evalInsLen = eval.getAttributeAsInt(GATKSVVCFConstants.INSERTED_SEQUENCE_LENGTH, 0);
                float diff = Math.abs(evalInsLen - truthInsLen);
                return 1 - insLenDiffScoreCoeff * Math.max(diff/evalInsLen, diff/truthInsLen);
            } else {
                return 0f;
            }
        } else if (evalType.equals(SimpleSVType.SupportedType.INV.name())) {

            if (evalType.equalsIgnoreCase(truthType)) { // NOTE: one can use SVLEN, but remember the truth set might report SVLEN as positive or 0 for INV
                SVInterval evalRange = new SVInterval(0, eval.getStart() - 1, eval.getEnd());
                SVInterval truthRange = new SVInterval(0, truth.getStart() - 1, truth.getEnd());
                float overlapLen = evalRange.overlapLen(truthRange);
                return Math.min(overlapLen/evalRange.getLength(), overlapLen/truthRange.getLength());
            } else {
                return 0f;
            }
        } else { // TODO: 6/4/18 all BND types will be considered discordant, for now
            return 0f;
        }
    }

    @Override
    public Object onTraversalSuccess() {
        printSummaryToStdErr();
        writeConcordanceToBedFile();
        return "SUCCESS";
    }

    // =================================================================================================================

    // TODO: 8/8/18 a nice feature to add is stratification by some features, e.g. genomic context, length, etc., and possibly even allowing users to specify the stratification key and levels
    private void printSummaryToStdErr() {

        // header
        System.err.print(String.format("%-25s:  ", "type"));
        for (final StructVarType type : StructVarType.values()) {System.err.print(String.format("%10s", type.name()));}
        System.err.println();

        final EnumMap<ConcordanceState, EnumMap<StructVarType, Long>> summaryByConcordanceThenBySvType =
                classifiedVariants.stream()
                        .collect(Collectors.groupingBy(EvalVariant::getConcordanceState,
                                () -> new EnumMap<>(ConcordanceState.class),
                                Collectors.groupingBy(classifiedEvalVariant ->
                                                StructVarType.valueOf(classifiedEvalVariant.getEval().getAttributeAsString(GATKSVVCFConstants.SVTYPE, "")),
                                        () -> new EnumMap<>(StructVarType.class),
                                        Collectors.counting())));
        summaryByConcordanceThenBySvType.forEach((concordanceState, countsBytype) -> {
            System.err.print(String.format("%-25s:  ", concordanceState.name()));
            final EnumMap<StructVarType, Long> zeroInitializedEnumMap = SVUtils.getZeroInitializedEnumMap(StructVarType.class);
            countsBytype.forEach(zeroInitializedEnumMap::put);
            zeroInitializedEnumMap.forEach((k, v) -> System.err.print(String.format("%10d", v)));
            System.err.print("\n");
        });

        // for FN values
        final EnumMap<StructVarType, Long> zeroInitializedEnumMap = SVUtils.getZeroInitializedEnumMap(StructVarType.class);
        System.err.print(String.format("%-25s:  ", ConcordanceState.FALSE_NEGATIVE.name()));
        Utils.stream(truthTree.iterator())
                .flatMap(e -> e.getValue().stream())
                .filter(classifiedTruth -> classifiedTruth.getConcordanceState().equals(ConcordanceState.FALSE_NEGATIVE))
                .collect(Collectors.groupingBy(variant -> StructVarType.valueOf(variant.getTruth().getAttributeAsString(GATKSVVCFConstants.SVTYPE, "")),
                        () -> new EnumMap<>(StructVarType.class), Collectors.counting()))
                .forEach(zeroInitializedEnumMap::put);
        zeroInitializedEnumMap.forEach((k,v) -> System.err.print(String.format("%10d", v)));

        System.err.println();
    }

    /**
     * BED file contains only: TP, FP, FILTERED_TRUE_NEGATIVE (record correctly filtered), and FILTERED_FALSE_NEGATIVE.
     *
     * False negatives are purposefully ignored because when first implemented, the truth set is from PacBio calls,
     * with a sensitivity that short read technology has no way of catching up to, generally speaking.
     */
    private void writeConcordanceToBedFile() {
        try (final OutputStreamWriter writer = new OutputStreamWriter(new BufferedOutputStream(
                BucketUtils.createFile(outPrefix + "concordance.bed")))) {
            for (final EvalVariant classifiedVariant : classifiedVariants) {

                final VariantContext eval = classifiedVariant.getEval();
                final String type = eval.getAttributeAsString(GATKSVVCFConstants.SVTYPE, "");

                if (type.equals(GATKSVVCFConstants.BREAKEND_STR))
                    continue;

                final String other;
                final float supportScore;
                if (classifiedVariant instanceof FalsePositive) {
                    other = ((FalsePositive) classifiedVariant).getReasonForFalsePositive().getAbbreviation();
                    supportScore = 0.0f; // TODO: 7/11/18 improve support score for FP, somtimes it's just border line support
                } else if (classifiedVariant instanceof FilteredTrueNegative) {
                    other = StringUtil.join(",", classifiedVariant.getEval().getFilters());
                    supportScore = 1.0f;
                } else {
//                        (classifiedVariant instanceof TruePositive || classifiedVariant instanceof FilteredFalseNegative)
                    other = ConcordanceState.TRUE_POSITIVE.getAbbreviation(); // to pad so that all fields have same number of fields
                    supportScore = ((BiVariant) classifiedVariant).getSupportScore();
                }
                final String concordance = classifiedVariant.getConcordanceState().getAbbreviation();
                final String evalID = eval.getID();
                final String contigNames = eval.hasAttribute(GATKSVVCFConstants.CONTIG_NAMES)
                        ? SVUtils.getAttributeAsStringStream(eval, GATKSVVCFConstants.CONTIG_NAMES).collect(Collectors.joining(","))
                        : "NA";
                writer.write(String.format("%s\t%d\t%d\t%s\t%.2f\n", eval.getContig(), eval.getStart(), eval.getEnd(),
                        concordance + ";" + type + ";" + other + ";" + evalID + ";" + contigNames, supportScore));
            }
        } catch (final IOException ioe) {
            throw new GATKException("Can't write BED file " + outPrefix + "concordance.bed", ioe);
        }
    }

    // TODO: 6/3/18 preferably, this should be in htsjdk
    private enum StructVarType {
        BND, INS, INV, DEL, DUP
    }
}
