package org.broadinstitute.hellbender.tools.walkers.sv;

import com.google.common.collect.Sets;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.StringUtil;
import htsjdk.tribble.bed.BEDFeature;
import htsjdk.tribble.bed.SimpleBEDFeature;
import htsjdk.variant.variantcontext.StructuralVariantType;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.commons.lang.StringUtils;
import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.engine.filters.VariantFilter;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SimpleSVType;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SvDiscoverFromLocalAssemblyContigAlignmentsSpark;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVIntervalTree;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVUtils;
import org.broadinstitute.hellbender.tools.walkers.validation.ConcordanceState;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import javax.annotation.Nonnull;
import javax.annotation.Nullable;
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
 *     <li>prints to screen a summary table containing summary on the TP/FP/FN/FTN/FFN numbers.</li>
 *     <li>writes a BED-4 file containing: TP, FP, FILTERED_TRUE_NEGATIVE (record correctly filtered), and FILTERED_FALSE_NEGATIVE</li>
 * </ul>
 *
 * <h3>Usage example</h3>
 * <pre>
 *   gatk EvaluateConcordanceOnNonGenotypedSVVCF \
 *     -T truth.vcf \
 *     -V evalCallSet.vcf \
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
public final class EvaluateConcordanceOnNonGenotypedSVVCF extends VariantWalker {

    public static final String CONFIDENCE_REGION_SHORT_NAME = "C";
    public static final String CONFIDENCE_REGION_LONG_NAME = "confidence";
    public static final String FEATURE_NAME_FOR_HIGH_CONFIDENCE = CONFIDENCE_REGION_LONG_NAME;

    public static final String MASK_REGION_SHORT_NAME = "M";
    public static final String MASK_REGION_LONG_NAME = "mask";
    public static final String FEATURE_NAME_FOR_MASKED_OUT = MASK_REGION_LONG_NAME;

    public static final String NON_CANONICAL_CHROMOSOMES_FILE_SHORT_NAME = "alt-tigs";
    public static final String NON_CANONICAL_CHROMOSOMES_FILE_LONG_NAME = "non-canonical-contig-names-file";

    public static final String EVALUATION_ON_NON_HUMAN_SHORT_NAME = "nonhs";
    public static final String EVALUATION_ON_NON_HUMAN_LONG_NAME = "non-human";

    public static final String PADDING_LONG_NAME = "PADDING";
    public static final String CONCORDANCE_RECIPROCAL_OVERLAP_LONG_NAME = "recip-ovp-frac";

    public static final String MAXIMUM_DIFFERENT_IN_INS_LEN_FOR_CONCORDANCE_LONG_NAME = "max-ins-len-diff-frac";
    public static final String MINIMUM_MAPPING_QUALITY_LONG_NAME = "min-map-qual";

    public static final String MINIMUM_ALIGNMENT_LENGTH_LONG_NAME = "min-aln-len";

    public static final String TRUTH_VARIANTS_SHORT_NAME = "T";
    public static final String TRUTH_VARIANTS_LONG_NAME = "truth";
    public static final String FEATURE_NAME_FOR_TRUTH = TRUTH_VARIANTS_LONG_NAME;

    public static final String CHR_Y_STRING_ARG_SHORT_NAME = "chry";
    public static final String CHR_Y_STRING_ARG_LONG_NAME = "chr-y-string";

    private static final Set<String> CHROMOSOME_Y_NAMES = Sets.newHashSet("Y", "chrY");

    private static final String NA_STRING = "NA";

    private static final String TEMPORARY_CONCORDANCE_STRING_FOR_TRUTH = "CONCORDANCE";

    public static final String NO_NEARBY_TRUTH = "NO_NEARBY_TRUTH";

    // ARGUMENTS =======================================================================================================

    /* file paths*/

    @Argument(doc = "path to a file holding high confidence region, if provided, variants out of this region will be skipped;" +
            " CURRENTLY UNUSED",
            shortName = CONFIDENCE_REGION_SHORT_NAME,
            fullName= CONFIDENCE_REGION_LONG_NAME,
            optional = true)
    private String highConfidence;

    @Argument(doc = "path to a file holding mask on region to be excluded, if provided, variants inside this region will be skipped",
            shortName = MASK_REGION_SHORT_NAME,
            fullName = MASK_REGION_LONG_NAME,
            optional = true)
    private String refRegionMaskFilePath;

    @Argument(doc = "file containing non-canonical chromosome names (e.g chrUn_KI270588v1) in the reference; " +
            "all chromosome will be considered canonical if not provided",
            shortName = NON_CANONICAL_CHROMOSOMES_FILE_SHORT_NAME,
            fullName = NON_CANONICAL_CHROMOSOMES_FILE_LONG_NAME,
            optional = true)
    private String nonCanonicalChromosomeNamesFile;

    @Argument(doc = "A VCF file containing truth variants",
            shortName = TRUTH_VARIANTS_SHORT_NAME,
            fullName = TRUTH_VARIANTS_LONG_NAME)
    private String truthVariantsFile;

    @Advanced
    @Argument(doc = "indicate evaluation is being performed on a non-human sample call set",
            shortName = EVALUATION_ON_NON_HUMAN_SHORT_NAME,
            fullName = EVALUATION_ON_NON_HUMAN_LONG_NAME,
            optional = true)
    private Boolean onNonHuman = false;

    @Advanced
    @Argument(doc = "string value for chromosome Y; takes value of \"Y\" or \"chrY\"; when left unprovided, " +
            "evaluation is skipped on chromosome Y; applicable only when \"non-human\" is false",
            shortName = CHR_Y_STRING_ARG_SHORT_NAME,
            fullName = CHR_Y_STRING_ARG_LONG_NAME,
            optional = true)
    private String chrYString;

    @Argument(doc = "output BED file holding information on evaluation call set variants",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME)
    private String outBed;

    /* Numeric value options*/

    @Argument(doc = "A fixed amount of padding for each VCF record from the eval call set. " +
            "Each individual record might be padded (more or less bases) as judged appropriate by the tool itself",
            fullName = PADDING_LONG_NAME,
            optional = true)
    private Integer fixedPadding = 50;

    @Argument(doc = "Minimum reciprocal overlap fraction for two concordant variants of the type DEL and INV",
            fullName = CONCORDANCE_RECIPROCAL_OVERLAP_LONG_NAME,
            optional = true)
    private Float recipOvpFracThreshold = 0.5f;

    @Argument(doc = "Maximum difference in length between two concordant variants or type INS and/or DUP, " +
            "i.e. the difference in their length cannot exceed this fraction of the shorter of the two",
            fullName = MAXIMUM_DIFFERENT_IN_INS_LEN_FOR_CONCORDANCE_LONG_NAME,
            optional = true)
    private Float insLenDiffFrac = 0.1f;

    @Advanced
    @Argument(doc = "Minimum mapping quality of the evidence contig mappings for a eval VCF record to be considered;" +
            " setting to zero turns off this filter",
            fullName = MINIMUM_MAPPING_QUALITY_LONG_NAME,
            optional = true)
    private Integer minMQ = 0;

    @Advanced
    @Argument(doc = "Minimum alignment length of the evidence contig mappings for a eval VCF record to be considered;" +
            " setting to zero turns off this filter",
            fullName = MINIMUM_ALIGNMENT_LENGTH_LONG_NAME,
            optional = true)
    private Integer minAlnLen = 0;

    // FIELDS ==========================================================================================================

    private SAMSequenceDictionary refSeqDict;

    private SVIntervalTree<String> maskedOut;
    private SVIntervalTree<String> highConfidenceRegions; // TODO: 9/4/18 this is not used yet in the current implementation, as we don't know what a high confidence region is for SV calls

    private SVIntervalTree<List<VariantContext>> truthTree;

    private Set<String> chromosomesOfInterest;

    // TODO: 6/4/18 this is initialized in a way that $ 1 - insLenDiffScoreCoeff * lenDiff/Math.min(evalLen, truthLen) < recipOvpFracThreshold $ for discordant insertion/duplication
    private Float insLenDiffScoreCoeff;

    private final List<EvalVariant> classifiedVariants = new ArrayList<>(50_000); // guess

    // INITIALIZATION ==================================================================================================

    @Override
    public final void onTraversalStart() {
        IOUtils.assertFileIsReadable(IOUtils.getPath(truthVariantsFile));

        // TODO: 8/31/18 add check that we are evaluating on the same sample

        refSeqDict = getBestAvailableSequenceDictionary();
        if (refRegionMaskFilePath != null) {
            maskedOut = buildIntervalTreeForMaskOrHighConfidenceRegion(refRegionMaskFilePath, refSeqDict, FEATURE_NAME_FOR_MASKED_OUT);
        } else {
            maskedOut = new SVIntervalTree<>();
        }

        if (highConfidence != null) {
            highConfidenceRegions = buildIntervalTreeForMaskOrHighConfidenceRegion(highConfidence, refSeqDict, FEATURE_NAME_FOR_HIGH_CONFIDENCE);
        } else {
            highConfidenceRegions = new SVIntervalTree<>();
        }

        buildTruth();

        insLenDiffScoreCoeff = (1f - recipOvpFracThreshold)/insLenDiffFrac;
        if (nonCanonicalChromosomeNamesFile != null) {
            chromosomesOfInterest = SVUtils.getCanonicalChromosomes(nonCanonicalChromosomeNamesFile, refSeqDict);
        } else {
            chromosomesOfInterest = refSeqDict.getSequences().stream().map(SAMSequenceRecord::getSequenceName)
                                              .collect(Collectors.toCollection(HashSet::new));
        }
        if (!onNonHuman && StringUtils.isEmpty(chrYString)) {
            chromosomesOfInterest.removeAll(CHROMOSOME_Y_NAMES);
        }
    }

    /**
     * Building an interval tree of high confidence regions from a BED file.
     * @param confidenceOrMaskRegionBedFilePath a file that follows the BED format specification; only the first 3 columns are used.
     * @param refSeqDict                        sequence dictionary for converting chromosome names to chromosome indexes
     * @param featureName                       feature name
     * @return                                  an interval tree of high confidence regions
     */
    private static SVIntervalTree<String> buildIntervalTreeForMaskOrHighConfidenceRegion(final String confidenceOrMaskRegionBedFilePath,
                                                                                         final SAMSequenceDictionary refSeqDict,
                                                                                         final String featureName) {
        final SVIntervalTree<String> result = new SVIntervalTree<>();
        try ( final FeatureDataSource<BEDFeature> regions =
                      new FeatureDataSource<>(confidenceOrMaskRegionBedFilePath, featureName, FEATURE_CACHE_LOOKAHEAD, SimpleBEDFeature.class) ) {
            regions.iterator().forEachRemaining(region ->
                result.put(new SVInterval(
                        refSeqDict.getSequenceIndex( region.getContig() ),
                        region.getStart(),
                        region.getEnd() + 1), "")
            );

            return result;
        }
    }

    /**
     * Note, certain coordinate conversions in this methods are confusing at first sight, because of VCF spec ambiguity.
     * First, for simple insertions and deletions, the spec REQUIREs the POS field to be the base before the deleted range
     * or the insertion site. But for duplications and inversions, the spec (version 4.3) reads (quote below)
     *
     * <p>
     *     For simple insertions and deletions in which either the REF or one of the ALT alleles would otherwise be null/empty,
     *     the REF and ALT Strings must include the base before the event (which must be reflected in the POS field),
     *     unless the event occurs at position 1 on the contig in which case it must include the base after the event;
     *     this padding base is not required (although it is permitted) for e.g.
     *     complex substitutions or other events where all alleles have at least one base represented in their Strings.
     *     If any of the ALT alleles is a symbolic allele (an angle-bracketed ID String “<ID>”)
     *     then the padding base is required and POS denotes the coordinate of the base preceding the polymorphism.
     * </p>
     *
     * All in all, SV call sets are not particularly spec-conformant, so we assume:
     * DUP and INV records, if available, are following the spec's requirement for symbolic allele variants.
     */
    private void buildTruth() {

        try ( final FeatureDataSource<VariantContext> truthVariants =
                new FeatureDataSource<>(truthVariantsFile, FEATURE_NAME_FOR_TRUTH, FEATURE_CACHE_LOOKAHEAD, VariantContext.class) ) {
            if (truthVariants.getSequenceDictionary() == null) {
                throw new UserException("Given truth variant VCF doesn't seem to have a valid reference sequence dictionary: " + truthVariantsFile);
            }

            truthTree = new SVIntervalTree<>();
            truthVariants.forEach(truth -> {

                if (truth.isFiltered()) return;

                final SVInterval truthInterval = convertSymbolicSVVCFRecord(truth, refSeqDict, 0, 0);
                final SVIntervalTree.Entry<List<VariantContext>> listEntry = truthTree.find(truthInterval);
                // pre-classify every truth call as FN, that is, without matching eval call, later could be modified accordingly by apply(...)
                final VariantContext truthDefaultsToFalseNegative = new VariantContextBuilder(truth)
                        .attribute(TEMPORARY_CONCORDANCE_STRING_FOR_TRUTH, ConcordanceState.FALSE_NEGATIVE.getAbbreviation()).make();
                final List<VariantContext> value;
                if (listEntry == null) {
                    value = Collections.singletonList(truthDefaultsToFalseNegative);
                } else {
                    value = new ArrayList<>(listEntry.getValue()); // listEntry is immutable
                    value.add(truthDefaultsToFalseNegative);
                }
                truthTree.put(truthInterval, value);
            });
        }
    }

    /**
     * For symbolic SV VCF records, the POS field has been shifted 1 base to the left, mandated by the spec
     * <p>
     *     For simple insertions and deletions in which either the REF or one of the ALT alleles would otherwise be null/empty,
     *     the REF and ALT Strings must include the base before the event (which must be reflected in the POS field),
     *     unless the event occurs at position 1 on the contig in which case it must include the base after the event;
     *     this padding base is not required (although it is permitted) for e.g.
     *     complex substitutions or other events where all alleles have at least one base represented in their Strings.
     *     If any of the ALT alleles is a symbolic allele (an angle-bracketed ID String “<ID>”)
     *     then the padding base is required and POS denotes the coordinate of the base preceding the polymorphism.
     * </p>
     */
    private static SVInterval convertSymbolicSVVCFRecord(@Nonnull final VariantContext symbolicSVCall,
                                                         @Nonnull final SAMSequenceDictionary sequenceDictionary,
                                                         final int frontPadding, final int rearPadding) {
        final int sequenceIndex = sequenceDictionary.getSequenceIndex(symbolicSVCall.getContig());
        if (sequenceIndex < 0)
            throw new IllegalArgumentException("Provided locatable: " + symbolicSVCall.toString() +
                    " doesn't seem to live on the chromosomes held in the provided dictionary: " + sequenceDictionary.toString());

        int start = Math.max(0, symbolicSVCall.getStart() + 1 - frontPadding);
        int end = Math.min(sequenceDictionary.getSequence(sequenceIndex).getSequenceLength(), symbolicSVCall.getEnd() + 1 + rearPadding);// + 1 is to count the fact that SVIntervals are right-open intervals
        return new SVInterval(sequenceIndex, start, end);
    }

    private static SimpleBEDFeature convertSymbolicSVVCFRecordToBED(@Nonnull final VariantContext symbolicSVCall) {
        return new SimpleBEDFeature(symbolicSVCall.getStart(), // remember, the POS has been shifted left by 1
                symbolicSVCall.getEnd(), // END has not been shifted, but BED format is right-open
                symbolicSVCall.getContig());
    }

    // TERMINATION =====================================================================================================

    @Override
    public Object onTraversalSuccess() {
        logSummary();
        writeConcordanceToBedFile();
        return "SUCCESS";
    }

    // TODO: 8/8/18 a nice feature to add is stratification by some features, e.g. genomic context, length, etc., and possibly even allowing users to specify the stratification key and levels
    private void logSummary() {

        logger.info("Evaluated a total number of " + classifiedVariants.size() + " variants");

        // header
        final StringBuilder headerPrint = new StringBuilder( String.format("%-25s:  ", "type") );
        for (final StructuralVariantType type : StructuralVariantType.values()) {headerPrint.append( String.format("%10s", type.name()) );}
        logger.info(headerPrint.toString());

        // collect summary by concordance state, then by variant type (note that FN are specifically skipped here)
        final EnumMap<ConcordanceState, EnumMap<StructuralVariantType, Long>> summaryByConcordanceThenBySvType =
                classifiedVariants.stream()
                        .collect(Collectors.groupingBy(EvalVariant::getAssignedConcordanceState,
                                () -> new EnumMap<>(ConcordanceState.class),
                                Collectors.groupingBy(classifiedEvalVariant ->
                                                StructuralVariantType.valueOf(classifiedEvalVariant.getVariant().getAttributeAsString(VCFConstants.SVTYPE, "")),
                                        () -> new EnumMap<>(StructuralVariantType.class),
                                        Collectors.counting())));
        summaryByConcordanceThenBySvType.forEach((concordanceState, countsBytype) -> {
            final StringBuilder concordanceStateAccountingLine = new StringBuilder(String.format("%-25s:  ", concordanceState.name()));
            final EnumMap<StructuralVariantType, Long> zeroInitializedEnumMap = SVUtils.getZeroInitializedEnumMap(StructuralVariantType.class);
            countsBytype.forEach(zeroInitializedEnumMap::put);
            zeroInitializedEnumMap.forEach((k, v) -> concordanceStateAccountingLine.append( String.format("%10d", v)) );
            logger.info(concordanceStateAccountingLine.toString());
        });

        // for FN values (motivation: using truth set from PacBio SV calls, the FN number is going to be high anyway for call sets based off short reads)
        final EnumMap<StructuralVariantType, Long> zeroInitializedEnumMap = SVUtils.getZeroInitializedEnumMap(StructuralVariantType.class);
        final StringBuilder falseNegativeLine = new StringBuilder(String.format("%-25s:  ", ConcordanceState.FALSE_NEGATIVE.name()));
        Utils.stream(truthTree.iterator())
                .flatMap(e -> e.getValue().stream())
                .filter(classifiedTruth -> classifiedTruth.getAttributeAsString(TEMPORARY_CONCORDANCE_STRING_FOR_TRUTH, "").equals(ConcordanceState.FALSE_NEGATIVE.getAbbreviation()))
                .collect(Collectors.groupingBy(truth -> StructuralVariantType.valueOf(truth.getAttributeAsString(VCFConstants.SVTYPE, "")),
                        () -> new EnumMap<>(StructuralVariantType.class), Collectors.counting()))
                .forEach(zeroInitializedEnumMap::put);
        zeroInitializedEnumMap.forEach((k,v) -> falseNegativeLine.append( String.format("%10d", v)) );

        logger.info(falseNegativeLine.toString());
    }

    /**
     * BED file contains only: TP, FP, FILTERED_TRUE_NEGATIVE (record correctly filtered), and FILTERED_FALSE_NEGATIVE.
     *
     * False negatives are purposefully ignored because when first implemented, the truth set is from PacBio calls,
     * with a sensitivity that short read technology has no way of catching up to, generally speaking.
     */
    private void writeConcordanceToBedFile() {
        try (final OutputStreamWriter writer = new OutputStreamWriter(new BufferedOutputStream(
                BucketUtils.createFile(outBed)))) {
            for (final EvalVariant classifiedVariant : classifiedVariants) {

                final VariantContext eval = classifiedVariant.getVariant();
                final String type = eval.getAttributeAsString(VCFConstants.SVTYPE, "");

                if (type.equals(GATKSVVCFConstants.BREAKEND_STR))
                    continue;

                final ConcordanceState assignedConcordanceState = classifiedVariant.getAssignedConcordanceState();
                final String id = eval.getID();

                final String other; // misc type of information
                final float supportScore;
                if ( assignedConcordanceState.equals(ConcordanceState.FALSE_POSITIVE) ) {
                    final NearbyTruth nearbyTruthMayBeNull = classifiedVariant.getNearbyTruthMayBeNull();
                    if (nearbyTruthMayBeNull == null) {
                        other = NO_NEARBY_TRUTH;
                        supportScore = 0;
                    } else {
                        other = nearbyTruthMayBeNull.getSupportType().name();
                        supportScore = nearbyTruthMayBeNull.supportScore;
                    }
                } else if ( assignedConcordanceState.equals(ConcordanceState.FILTERED_TRUE_NEGATIVE) ) {
                    other = StringUtil.join(",", classifiedVariant.getVariant().getFilters());
                    supportScore = 1.0f;
                } else if ( assignedConcordanceState.equals(ConcordanceState.FILTERED_FALSE_NEGATIVE) ){
                    other = StringUtil.join(",", classifiedVariant.getVariant().getFilters());
                    supportScore = classifiedVariant.getNearbyTruthMayBeNull().getSupportScore();
                } else {
                    other = ConcordanceState.TRUE_POSITIVE.getAbbreviation(); // to pad so that all fields have same number of fields
                    supportScore = classifiedVariant.getNearbyTruthMayBeNull().getSupportScore();
                }

                final String contigNames = eval.hasAttribute(GATKSVVCFConstants.CONTIG_NAMES)
                        ? SVUtils.getAttributeAsStringStream(eval, GATKSVVCFConstants.CONTIG_NAMES).collect(Collectors.joining(","))
                        : NA_STRING;
                final SimpleBEDFeature simpleBEDFeature = convertSymbolicSVVCFRecordToBED(eval);
                writer.write(String.format("%s\t%d\t%d\t%s\t%.2f\n", simpleBEDFeature.getContig(), simpleBEDFeature.getStart(), simpleBEDFeature.getEnd(),
                        assignedConcordanceState.getAbbreviation() + ";" + type + ";" + other + ";" + id + ";" + contigNames, supportScore));
            }
        } catch (final IOException ioe) {
            throw new GATKException("Can't write BED file " + outBed, ioe);
        }
    }

    // TRAVERSAL =======================================================================================================

    /**
     * A 2-step evaluation mechanism:
     * <ul>
     *     <li> first test if {@code eval} (padded by {@link #getPaddedSvInterval(VariantContext)}) overlaps with any record in the truth tree </li>
     *     <li> further test concordance by {@link #getNearbyTruth(VariantContext)} assuming there's such overlap </li>
     * </ul>
     */
    @Override
    public final void apply(final VariantContext eval, final ReadsContext readsContext,
                            final ReferenceContext refContext, final FeatureContext featureContext) {

        final boolean filtered = eval.isFiltered();

        final NearbyTruth nearbyTruth = getNearbyTruth(eval);
        final ConcordanceState assignedConcordanceState;
        if ( nearbyTruth != null && nearbyTruth.getSupportType().equals(NearbyTruth.SupportType.GOOD_SUPPORT) ){
            assignedConcordanceState = filtered ? ConcordanceState.FILTERED_FALSE_NEGATIVE : ConcordanceState.TRUE_POSITIVE;
        } else {
            assignedConcordanceState = filtered ? ConcordanceState.FILTERED_TRUE_NEGATIVE : ConcordanceState.FALSE_POSITIVE;
        }
        classifiedVariants.add( new EvalVariant(eval, nearbyTruth, assignedConcordanceState) );
    }

    /**
     * Pad provided eval call accordingly for overlapping with the truth interval tree.
     */
    private SVInterval getPaddedSvInterval(final VariantContext eval) {

        final int homologyLength = eval.getAttributeAsInt(GATKSVVCFConstants.HOMOLOGY_LENGTH, 0);

        int insTypeRightPadding = 0;
        final String svType = eval.getAttributeAsString(VCFConstants.SVTYPE, "");
        if (svType.equals(GATKSVVCFConstants.SYMB_ALT_ALLELE_INS) || svType.equals(GATKSVVCFConstants.SYMB_ALT_ALLELE_DUP))
            insTypeRightPadding = eval.getAttributeAsInt(GATKSVVCFConstants.SVLEN, 0);

        // we only pad to the right with homology because we expect variants to follow left-align convention (GATK-SV does this)
        final int leftPadding = fixedPadding;
        final int rightPadding = fixedPadding + homologyLength + insTypeRightPadding;

        return convertSymbolicSVVCFRecord(eval, refSeqDict, leftPadding, rightPadding);
    }

    // UPFRONT FILTERS APPLIED TO CALLS

    @Override
    protected final VariantFilter makeVariantFilter() {

        VariantFilter finalFilter = filterOutVariantOnUninterestingChromosomes(chromosomesOfInterest)
                                    .and(filterOutSmallVariant());
        if (minMQ != 0)
            finalFilter = finalFilter.and(new SvDiscoverFromLocalAssemblyContigAlignmentsSpark.SVMappingQualityFilter(minMQ));
        if (minAlnLen != 0)
            finalFilter = finalFilter.and(new SvDiscoverFromLocalAssemblyContigAlignmentsSpark.SVAlignmentLengthFilter(minAlnLen));
        if (maskedOut.size() != 0)
            finalFilter = finalFilter.and(filterOutVariantContainedInMaskedOutSites(maskedOut, refSeqDict));
        return finalFilter;
    }

    private static VariantFilter filterOutVariantOnUninterestingChromosomes(final Set<String> chromosomesOfInterest) {
        return variantContext -> chromosomesOfInterest.contains(variantContext.getContig());
    }

    /**
     * Filters out variants that are completely contained in masked out sites,
     * that is, simply overlap will NOT necessarily filter out the variant.
     */
    private static VariantFilter filterOutVariantContainedInMaskedOutSites(final SVIntervalTree<? extends Object> maskedOutSites,
                                                                           final SAMSequenceDictionary refSeqDict) {
        return variantContext -> {
            final SVInterval svInterval = convertSymbolicSVVCFRecord(variantContext, refSeqDict, 0, 0);
            if ( maskedOutSites.hasOverlapper(svInterval) ) {
                Iterator<? extends SVIntervalTree.Entry<?>> overlappers = maskedOutSites.overlappers(svInterval);
                while (overlappers.hasNext()) {
                    if (overlappers.next().getInterval().contains(svInterval))
                        return false;
                }
                return true;
            } else
                return true;
        };
    }

    /**
     * Occasionally, callers emit small variants but we don't want to look at them.
     */
    private static VariantFilter filterOutSmallVariant() {
        return variantContext -> {
            final int svLen = Math.abs(variantContext.getAttributeAsInt(GATKSVVCFConstants.SVLEN, 0));
            final String svType = variantContext.getAttributeAsString(VCFConstants.SVTYPE, "");
            if (svType.equals(SimpleSVType.SupportedType.DEL.name())) { // some tiny deletion calls are actually RPL calls with micro deletions
                if ( (!variantContext.hasAttribute(GATKSVVCFConstants.INSERTED_SEQUENCE)) && svLen < StructuralVariationDiscoveryArgumentCollection.STRUCTURAL_VARIANT_SIZE_LOWER_BOUND) { // no ins_seq AND tiny
                    return false;
                }
            } else if (svType.equals(SimpleSVType.SupportedType.INS.name()) || svType.equals(SimpleSVType.SupportedType.DUP.name())) {
                if (svLen < StructuralVariationDiscoveryArgumentCollection.STRUCTURAL_VARIANT_SIZE_LOWER_BOUND) {
                    return false;
                }
            }
            return true;
        };
    }

    // CONCORDANCE

    private static final class EvalVariant {
        private final VariantContext variant;
        private final NearbyTruth nearbyTruthMayBeNull;
        private ConcordanceState assignedConcordanceState;

        // when {@code nearbyTruthMayBeNull} is null, it means no truth variant is near by
        EvalVariant(@Nonnull final VariantContext variant, @Nullable final NearbyTruth nearbyTruthMayBeNull,
                    @Nonnull final ConcordanceState assignedConcordanceState) {
            this.variant = variant;
            this.nearbyTruthMayBeNull = nearbyTruthMayBeNull;
            this.assignedConcordanceState = assignedConcordanceState;
        }

        private VariantContext getVariant() {
            return variant;
        }

        private NearbyTruth getNearbyTruthMayBeNull() {
            return nearbyTruthMayBeNull;
        }

        private ConcordanceState getAssignedConcordanceState() {
            return assignedConcordanceState;
        }
    }

    /**
     * A {@link EvalVariant} can have several nearby truth variants,
     * (some will have none, meaning no truth variant is in the vicinity of it)
     * each might offer a different kind of support,
     * some might be of incompatible type (though nearby),
     * some might be of compatible type but offers low support
     *  (e.g. del/inv of very different range, ins/dup of very different inserted seq).
     * This class is to symbolize that.
     */
    private static final class NearbyTruth {

        enum SupportType {
            INCOMPATIBLE_TYPE,
            LOW_SUPPORT,
            GOOD_SUPPORT
        }

        private final VariantContext truthVar;
        private final float supportScore;
        private final SupportType supportType;

        NearbyTruth(@Nonnull final VariantContext truthVar, final float supportScore, @Nonnull final SupportType supportType) {
            this.truthVar = truthVar;
            this.supportScore = supportScore;
            this.supportType = supportType;
        }

        VariantContext getTruthVar() {
            return truthVar;
        }

        float getSupportScore() {
            return supportScore;
        }

        SupportType getSupportType() {
            return supportType;
        }
    }

    /**
     * Given an eval variant, return concordance truth variant, if any.
     * @param eval  variant to be evaluated
     * @return      {@code null} if no near by truth exists
     */
    private NearbyTruth getNearbyTruth(final VariantContext eval) {

        NearbyTruth bestNearbyTruth = null;

        final Iterator<SVIntervalTree.Entry<List<VariantContext>>> nearbyTruths = truthTree.overlappers(getPaddedSvInterval(eval));

        float maxConcordanceScore = 0f;
        while (nearbyTruths.hasNext()) {
            final List<VariantContext> truthVariants = new ArrayList<>( nearbyTruths.next().getValue() );
            for (int i = 0; i < truthVariants.size(); ++i ) {

                final VariantContext truthVarCtx = truthVariants.get(i);
                if ( ! typesCompatible(eval, truthVarCtx) && bestNearbyTruth == null) {
                    bestNearbyTruth = new NearbyTruth(truthVarCtx, 0, NearbyTruth.SupportType.INCOMPATIBLE_TYPE); // TODO: 9/6/18 here we have some ambiguity: if none of the nearby truths are of compatible type, this always sets the incompatible nearby truth as the one of the highest coordinate, but we are not currently using this truth downstream, so it's ok for now
                } else { // at least types are compatible
                    final float concordanceScore = computeConcordanceScore(eval, truthVarCtx, null);

                    // TODO: 6/3/18 test if heuristic overlap-based concordance score is over a hard threshold, in the long run, haplotype concordance should be evaluated instead
                    final boolean concordant = recipOvpFracThreshold <= concordanceScore;

                    if (concordant) {
                        if (concordanceScore > maxConcordanceScore) { // conditionally update to best support
                            final VariantContextBuilder variantContextBuilder = new VariantContextBuilder(truthVarCtx).rmAttribute(TEMPORARY_CONCORDANCE_STRING_FOR_TRUTH);
                            final VariantContext updatedTruthVariant = eval.isFiltered() ? variantContextBuilder.attribute(TEMPORARY_CONCORDANCE_STRING_FOR_TRUTH, ConcordanceState.FILTERED_FALSE_NEGATIVE.getAbbreviation()).make()
                                                                                         : variantContextBuilder.attribute(TEMPORARY_CONCORDANCE_STRING_FOR_TRUTH, ConcordanceState.TRUE_POSITIVE.getAbbreviation()).make();
                            truthVariants.set(i, updatedTruthVariant);

                            bestNearbyTruth = new NearbyTruth(truthVarCtx, concordanceScore, NearbyTruth.SupportType.GOOD_SUPPORT);
                            maxConcordanceScore = concordanceScore;
                        }
                    } else if (bestNearbyTruth == null || bestNearbyTruth.getSupportType() != NearbyTruth.SupportType.GOOD_SUPPORT){
                        bestNearbyTruth = new NearbyTruth(truthVarCtx, concordanceScore, NearbyTruth.SupportType.LOW_SUPPORT);
                    }
                }
            }
        }
        return bestNearbyTruth;
    }

    private boolean typesCompatible(final VariantContext eval, final VariantContext truth) {
        final StructuralVariantType evalType = StructuralVariantType.valueOf(eval.getAttributeAsString(VCFConstants.SVTYPE, ""));
        final StructuralVariantType truthType = StructuralVariantType.valueOf(truth.getAttributeAsString(VCFConstants.SVTYPE, ""));

        if ( evalType.equals(truthType) ) {
            return true;
        } else if (evalType.equals(StructuralVariantType.CNV)){
            return truthType.equals(StructuralVariantType.DEL) || truthType.equals(StructuralVariantType.DUP) || truthType.equals(StructuralVariantType.INS);
        } else if (truthType.equals(StructuralVariantType.CNV)){
            return evalType.equals(StructuralVariantType.DEL) || evalType.equals(StructuralVariantType.DUP) || evalType.equals(StructuralVariantType.INS);
        } else {
            if (evalType.equals(StructuralVariantType.DUP) && truthType.equals(StructuralVariantType.INS)) {
                return true;
            }
            return evalType.equals(StructuralVariantType.INS) && truthType.equals(StructuralVariantType.DUP);
        }
    }

    /**
     * Computes a score of concordance between the eval and truth.
     * This is particularly useful for variants that spans a large interval, i.e. CNV and SV calls.
     * >= 0.5 is considered concordant.
     */
    private float computeConcordanceScore(final VariantContext eval, final VariantContext truth, final ReferenceContext context) {
        final StructuralVariantType evalType = StructuralVariantType.valueOf( eval.getAttributeAsString(VCFConstants.SVTYPE, "") );
        final StructuralVariantType truthType = StructuralVariantType.valueOf( truth.getAttributeAsString(VCFConstants.SVTYPE, "") );

        // TODO: 6/3/18 we only consider insertion length but DO NOT YET MATCH inserted sequence (we should, in the long run, but that involves haplotype reconstruction)
        if (evalType.equals(StructuralVariantType.INS) || evalType.equals(StructuralVariantType.DUP)) {
            if (truthType.equals(StructuralVariantType.INS)) { // TODO: 6/3/18 truth type doesn't call DUP, instead type them as INS
                int evalLen = eval.getAttributeAsInt(GATKSVVCFConstants.SVLEN, 0);
                int truthLen = truth.getAttributeAsInt(GATKSVVCFConstants.SVLEN, 0);
                float diff = Math.abs(evalLen - truthLen);
                return Math.max( 0f, 1 - insLenDiffScoreCoeff * Math.max(diff/evalLen, diff/truthLen) ); // cap score from below at 0.0; negative value can happen when the truth reported SVLEN is much higher than eval reported SVLEN
            } else
                return 0f;
        } else if (evalType.equals(StructuralVariantType.DEL)){

            if (truthType.equals(StructuralVariantType.DEL)) { // NOTE: one can use SVLEN, but remember the truth set might report positive SVLEN for DEL
                final SVInterval truthRange = convertSymbolicSVVCFRecord(truth, refSeqDict, 0, 0);
                final SVInterval evalRange = convertSymbolicSVVCFRecord(eval, refSeqDict, 0, 0);
                float overlapLen = evalRange.overlapLen(truthRange);
                return Math.min(overlapLen/evalRange.getLength(), overlapLen/truthRange.getLength());
            } else if (truthType.equals(StructuralVariantType.INS) && eval.hasAttribute(GATKSVVCFConstants.INSERTED_SEQUENCE_LENGTH)) { // some call set might report replacement as deletion with insertion
                final int truthInsLen = truth.getAttributeAsInt(GATKSVVCFConstants.SVLEN, 0);
                final int evalInsLen = eval.getAttributeAsInt(GATKSVVCFConstants.INSERTED_SEQUENCE_LENGTH, 0);
                float diff = Math.abs(evalInsLen - truthInsLen);
                return 1 - insLenDiffScoreCoeff * Math.max(diff/evalInsLen, diff/truthInsLen);
            } else {
                return 0f;
            }
        } else if (evalType.equals(StructuralVariantType.INV)) {

            if (evalType.equals(truthType)) { // NOTE: one can use SVLEN, but remember the truth set might report SVLEN as positive or 0 for INV
                final SVInterval truthRange = convertSymbolicSVVCFRecord(truth, refSeqDict, 0, 0);
                final SVInterval evalRange = convertSymbolicSVVCFRecord(eval, refSeqDict, 0, 0);
                float overlapLen = evalRange.overlapLen(truthRange);
                return Math.min(overlapLen/evalRange.getLength(), overlapLen/truthRange.getLength());
            } else {
                return 0f;
            }
        } else { // TODO: 6/4/18 all BND and CNV types will be considered discordant, for now
            return 0f;
        }
    }
}
