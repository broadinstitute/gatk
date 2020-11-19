package org.broadinstitute.hellbender.tools.walkers.variantutils;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.argumentcollections.DbsnpArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.gnarlyGenotyper.GnarlyGenotyper;
import org.broadinstitute.hellbender.utils.logging.OneShotLogger;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import org.broadinstitute.hellbender.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.hellbender.utils.reference.ReferenceUtils;
import picard.cmdline.programgroups.VariantEvaluationProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.VariantWalker;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.*;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.*;


/**
 * Validate a VCF file with a strict set of criteria
 *
 * <p> This tool is designed to validate the adherence of a file to VCF format. The tool will validate .g.vcf GVCF
 * format files as well. For VCF specifications, see
 * <a href='https://samtools.github.io/hts-specs/'>https://samtools.github.io/hts-specs/</a>.
 * Besides standard adherence to the VCF specification, this tool performs additional strict validations to ensure
 * that the information contained within the file is correctly encoded. These include:
 * </p>
 *
 * <ul>
 *   <li><b>REF</b> - correctness of the reference base(s)</li>
 *   <li><b>CHR_COUNTS</b> - accuracy of AC and AN values</li>
 *   <li><b>IDS</b> - tests against rsIDs when a dbSNP file is provided (requires a dbsnp VCF provided via `--dbsnp`).</li>
 *   <li><b>ALLELES</b> - that all alternate alleles are present in at least one sample</li>
 * </ul>
 *
 * <p>
 *     By default the tool applies all the strict validations unless you indicate which one should be
 *     excluded using `--validation-type-to-exclude`. You can exclude as many types as you want. Furthermore, you
 *     can exclude all strict validations with the special code `ALL`. In this case the tool will only test for
 *     adherence to the VCF specification.
 * </p>
 *
 * <h3>Input</h3>
 * <p>
 * A VCF file to validate.
 * </p>
 *
 * <h3>Usage examples</h3>
 *
 * <h4>Minimally validate a file for adherence to VCF format:</h4>
 * gatk ValidateVariants \
 *     -V cohort.vcf.gz
 *
 * <h4>Validate a GVCF for adherence to VCF format, including REF allele match:</h4>
 * gatk ValidateVariants \
 *     -V sample.g.vcf.gz \
 *     -R reference.fasta
 *     -gvcf
 *
 * <h4>To perform VCF format and all strict validations: </h4>
 * <pre>
 * gatk ValidateVariants \
 *   -R ref.fasta \
 *   -V input.vcf \
 *   --dbsnp dbsnp.vcf
 * </pre>
 *
 * <h4>To perform only VCF format tests:</h4>
 * <pre>
 * gatk ValidateVariants
 *   -R ref.fasta \
 *   -V input.vcf \
 *   --validation-type-to-exclude ALL
 * </pre>
 *
 * <h4>To perform all validations except the strict `ALLELE` validation:</h4>
 * <pre>
 * gatk ValidateVariants \
 *   -R ref.fasta \
 *   -V input.vcf \
 *   --validation-type-to-exclude ALLELES \
 *   --dbsnp dbsnp.vcf
 * </pre>
 *
 */
@CommandLineProgramProperties(
        summary = "Validates a VCF file with an extra strict set of criteria.",
        oneLineSummary = "Validate VCF",
        programGroup = VariantEvaluationProgramGroup.class
)
@DocumentedFeature
public final class ValidateVariants extends VariantWalker {
    static final Logger logger = LogManager.getLogger(ValidateVariants.class);
    static final OneShotLogger oneShotLogger = new OneShotLogger(ValidateVariants.class);

    public static final String GVCF_VALIDATE = "validate-GVCF";
    public static final String DO_NOT_VALIDATE_FILTERED_RECORDS = "do-not-validate-filtered-records";

    public enum ValidationType {

        /**
         * Check whether the reported reference base in the VCF is the same as the corresponding base in the
         * actual reference.
         */
        REF,

        /**
         * Checks whether the variant IDs exists, only relevant if the user indicates a DBSNP vcf file (see {@link #dbsnp}).
         */
        IDS,

        /**
         * Check whether all alternative alleles participate in a genotype call of at least on sample.
         */
        ALLELES,

        /**
         * Check that the AN and AC annotations are consistent with the number of calls, alleles and then number these
         * are called across samples.
         */
        CHR_COUNTS,

        /**
         * Check that each genotype has a GT and AD and (for sites with no more than 6 alt alleles) PLs and GQ
         */
        GNARLY_INPUT,

        /**
         * Check that each variant has required VQSR annotations, including rank sums if heterozygous genotypes are present
         */
        VQSR_INPUT,

        /**
         * Check that allele-specific annotations have the right number of values based on the listed ALTs
         */
        AS_ANNOTATIONS,

        /**
         * Includes reference checking (if reference is supplied), dbSNP (if dbSNP is supplied), alt alleles and chromosome counts
         */
        ALL;

        /**
         * Unmodifiable set of concrete validation types.
         *
         * <p>These are all types except {@link #ALL}.</p>
         */
        public static final Set<ValidationType> CONCRETE_TYPES;


        public static final Set<ValidationType> DEFAULT_TYPES;

        static {
            final Set<ValidationType> cts = new LinkedHashSet<>(values().length - 1);
            for (final ValidationType v : values()) {
                if (v != ALL)
                    cts.add(v);
            }
            CONCRETE_TYPES = Collections.unmodifiableSet(cts);
        }

        static {
            DEFAULT_TYPES = new LinkedHashSet<>(CONCRETE_TYPES);
            DEFAULT_TYPES.remove(ValidationType.GNARLY_INPUT);
            DEFAULT_TYPES.remove(ValidationType.VQSR_INPUT);
            DEFAULT_TYPES.remove(ValidationType.AS_ANNOTATIONS);
        }
    }

    /**
     * Contains final set of validation to apply.
     */
    boolean[] validationsToPerform;

    @ArgumentCollection
    DbsnpArgumentCollection dbsnp = new DbsnpArgumentCollection();

    @Argument(fullName = "validation-type-to-exclude",
            shortName = "Xtype",
            doc = "which validation type to exclude from a full strict validation",
            optional = true)
    List<ValidationType> excludeTypes = new ArrayList<>();

    @Argument(fullName = "validation-type-to-include",
            shortName = "optional-type",
            doc = "which optional validation type to include in validation (e.g. VQSR_INPUT, AS_ANNOTATIONS)",
            optional = true)
    List<ValidationType> includeTypes = new ArrayList<>();

    /**
     * By default, even filtered records are validated.
     */
    @Argument(fullName = DO_NOT_VALIDATE_FILTERED_RECORDS,
            shortName = "do-not-validate-filtered-records",
            doc = "skip validation on filtered records",
            optional = true,
            mutex = GVCF_VALIDATE)
    Boolean DO_NOT_VALIDATE_FILTERED = false;

    @Argument(fullName = "warn-on-errors",
            shortName = "warn-on-errors",
            doc = "just emit warnings on errors instead of terminating the run at the first instance",
            optional = true)
    Boolean WARN_ON_ERROR = false;

    /**
     *  Validate this file as a gvcf. In particular, every variant must include a <NON_REF> allele, and that
     *  every base in the territory under consideration is covered by a variant (or a reference block).
     *  If you specifed intervals (using -L or -XL) to restrict analysis to a subset of genomic regions,
     *  those intervals will need to be covered in a valid gvcf.
     */
    @Argument(fullName = GVCF_VALIDATE,
            shortName = "gvcf",
            doc = "Validate this file as a GVCF",
            optional = true,
            mutex = DO_NOT_VALIDATE_FILTERED_RECORDS)
    Boolean VALIDATE_GVCF = false;

    /**
     * Fail if GVCF contains positions that overlap. Overlapping variants are allowed because they are assumed to be on
     * separate haplotypes.
     */
    @Argument(fullName = "fail-gvcf-on-overlap",
            shortName = "no-overlaps",
            doc = "Fail if GVCF contains positions that overlap (except two variants overlapping)",
            optional = true,
            mutex = DO_NOT_VALIDATE_FILTERED_RECORDS)
    Boolean FAIL_ON_OVERLAP = false;

    /**
     * If using VQSR validation, validate the subset of annotations used for exomes
     */
    @Argument(fullName = "exome-input",
            shortName = "exome",
            doc = "Validate this file against expected annotations for VQSR for exomes",
            optional = true)
    Boolean EXOME_INPUT = false;

    @Argument(fullName = "max-alt-alleles",
            shortName = "max-alt-alleles",
            doc = "Maximum number of alt alleles for which PLs should be reported; used in GNARLY_INPUT validation mode",
            optional = true)
    int MAX_ALT_ALLELES = 5;

    private GenomeLocSortedSet genomeLocSortedSet;
    private CachingIndexedFastaSequenceFile referenceReader;

    // information to keep track of when validating a GVCF
    private SimpleInterval previousInterval;
    private int previousStart = -1;
    private String previousContig = null;
    private boolean previousIntervalIsReference = true;
    private boolean sawOverlap = false;
    private SimpleInterval firstOverlap;

    private static final List<String> requiredVQSRAnnotationKeys = Arrays.asList("MQ", "QD", "MQRankSum", "ReadPosRankSum", "FS", "SOR");  //DP is not required for exomes
    private static final List<String> requiredRawVQSRAnnotationKeys = Arrays.asList("RAW_MQandDP", "MQRankSum", "ReadPosRankSum");  //DP is not required for exomes
    private static final List<String> requiredAlleleSpecificVQSRAnnotationKeys = Arrays.asList("AS_MQ", "AS_QD", "AS_MQRankSum", "AS_ReadPosRankSum", "AS_FS", "AS_SOR");  //DP is not required for exomes, but should appear anyway
    private static final List<String> requiredRawASVQSRAnnotationKeys = Arrays.asList(GATKVCFConstants.AS_RAW_RMS_MAPPING_QUALITY_KEY, GATKVCFConstants.AS_RAW_MAP_QUAL_RANK_SUM_KEY, GATKVCFConstants.AS_RAW_READ_POS_RANK_SUM_KEY,
                                                                            GATKVCFConstants.AS_SB_TABLE_KEY);

    @Override
    public void onTraversalStart() {
        if (VALIDATE_GVCF) {
            final SAMSequenceDictionary seqDictionary = getBestAvailableSequenceDictionary();

            if (seqDictionary == null)
                throw new UserException("Validating a GVCF requires a sequence dictionary but no dictionary was able to be constructed from your input.");

            genomeLocSortedSet = new GenomeLocSortedSet(new GenomeLocParser(seqDictionary));

            if (hasReference()) {
                final SAMSequenceDictionary inputDict = getSequenceDictionaryForDrivingVariants();
                final SAMSequenceDictionary refDict = getReferenceDictionary();
                if (inputDict.isSameDictionary(refDict)) {
                    logger.warn("GVCF sequence dictionary does not match the sequence dictionary for the supplied reference.");
                }
            }
        }
        validationsToPerform = calculateValidationTypesToApply(excludeTypes);


        //warn user if certain requested validations cannot be done due to lack of arguments
        if(dbsnp.dbsnp == null && validationsToPerform[ValidationType.IDS.ordinal()])
        {
            logger.warn("IDS validation cannot be done because no DBSNP file was provided");
            logger.warn("Other possible validations will still be performed");
        }
        if(!hasReference() && validationsToPerform[ValidationType.REF.ordinal()])
        {
            logger.warn("REF validation cannot be done because no reference file was provided");
            logger.warn("Other possible validations will still be performed");
        }
        if(hasReference()) {
            referenceReader = ReferenceUtils.createReferenceReader(referenceArguments.getReferenceSpecifier());
        }
    }

    @Override
    public void apply(final VariantContext vc, final ReadsContext readsContext, final ReferenceContext ref, final FeatureContext featureContext) {
        if (DO_NOT_VALIDATE_FILTERED && vc.isFiltered()) {
            return;
        }
        // get the true reference allele
        final Allele reportedRefAllele = vc.getReference();
        final int refLength = reportedRefAllele.length();

        final Allele observedRefAllele = hasReference() ? Allele.create(ReferenceUtils.getRefBasesAtPosition(referenceReader, vc.getContig(), vc.getStart(), refLength)) : null;

        final Set<String> rsIDs = getRSIDs(featureContext);

        if (VALIDATE_GVCF) {
            final SimpleInterval refInterval = ref.getInterval();

            validateVariantsOrder(vc);

            final boolean thisIntervalIsReference =  vc.getAlternateAlleles().size() == 1 && vc.getAlternateAllele(0).equals(Allele.NON_REF_ALLELE);

            // GenomeLocSortedSet will automatically merge intervals that are overlapping when setting `mergeIfIntervalOverlaps`
            // to true.  In a GVCF most blocks are adjacent to each other so they wouldn't normally get merged.  We check
            // if the current record is adjacent to the previous record and "overlap" them if they are so our set is as
            // small as possible while still containing the same bases.
            if (previousInterval != null && previousInterval.overlapsWithMargin(refInterval, 0) &&
                    (previousIntervalIsReference || thisIntervalIsReference)) {
                logger.warn("Current interval " + refInterval.toString() + " overlaps previous interval ending at " + previousInterval.getEnd());
                sawOverlap = true;
                firstOverlap = refInterval;
            }
            final int start = (previousInterval != null && previousInterval.overlapsWithMargin(refInterval, 1)) ?
                    previousInterval.getStart() : refInterval.getStart();
            final int end = (previousInterval != null && previousInterval.overlapsWithMargin(refInterval, 1)) ?
                    Math.max(previousInterval.getEnd(), vc.getEnd()) : vc.getEnd();
            final GenomeLoc possiblyMergedGenomeLoc = genomeLocSortedSet.getGenomeLocParser().createGenomeLoc(refInterval.getContig(), start, end);
            genomeLocSortedSet.add(possiblyMergedGenomeLoc, true);

            previousInterval = new SimpleInterval(possiblyMergedGenomeLoc);
            previousStart = vc.getStart();
            previousIntervalIsReference = thisIntervalIsReference;
            validateGVCFVariant(vc);
        }
        applyValidations(vc, reportedRefAllele, observedRefAllele, rsIDs);
    }

    @Override
    public Object onTraversalSuccess() {
        if (VALIDATE_GVCF) {
            final GenomeLocSortedSet intervalArgumentGenomeLocSortedSet;
            final SAMSequenceDictionary seqDictionary = getBestAvailableSequenceDictionary();

            if (intervalArgumentCollection.intervalsSpecified()){
                intervalArgumentGenomeLocSortedSet = GenomeLocSortedSet.createSetFromList(genomeLocSortedSet.getGenomeLocParser(), IntervalUtils.genomeLocsFromLocatables(genomeLocSortedSet.getGenomeLocParser(), intervalArgumentCollection.getIntervals(seqDictionary)));
            } else {
                intervalArgumentGenomeLocSortedSet = GenomeLocSortedSet.createSetFromSequenceDictionary(seqDictionary);
            }

            final GenomeLocSortedSet uncoveredIntervals = intervalArgumentGenomeLocSortedSet.subtractRegions(genomeLocSortedSet);
            if (uncoveredIntervals.coveredSize() > 0) {
                final UserException e = new UserException.ValidationFailure("A GVCF must cover the entire region. Found " + uncoveredIntervals.coveredSize() +
                        " loci with no VariantContext covering it. The first uncovered segment is:" +
                        uncoveredIntervals.iterator().next());
                throwOrWarn(e);
            }
            if (FAIL_ON_OVERLAP && sawOverlap) {
                final UserException e = new UserException.ValidationFailure("This GVCF contained overlapping reference blocks.  The first overlapping interval is " +
                        firstOverlap.toString());
                throwOrWarn(e);
            }
        }
        return null;
    }

    /*
     *  Returns the list of RSIDs overlapping the current variant that we're walking over.
     *  If there's no RSID or if there was not dbsnp file passed in as an argument,
     *  an empty set is returned (and then no validation is performed, see applyValidations.
     */
    private Set<String> getRSIDs(FeatureContext featureContext) {
        Set<String> rsIDs = new LinkedHashSet<>();
        for (VariantContext rsID : featureContext.getValues(dbsnp.dbsnp)) {
            rsIDs.addAll(Arrays.asList(rsID.getID().split(VCFConstants.ID_FIELD_SEPARATOR)));
        }
        return rsIDs;
    }

    /**
     * Given the validation type and exclusion type, calculate the final set of type to validate.
     * @param excludeTypes types to exclude.
     *
     * @return the final set of type to validate. May be empty.
     */
    private boolean[] calculateValidationTypesToApply(final List<ValidationType> excludeTypes) {

        //creates local, temp list so that original list provided by user doesn't get modified
        final List<ValidationType> excludeTypesTemp = new ArrayList<>(excludeTypes);
        final Set<ValidationType> includeTypes;
        if (VALIDATE_GVCF && !excludeTypesTemp.contains(ValidationType.ALLELES)) {
            // Note: in a future version allele validation might be OK for GVCFs, if that happens
            // this will be more complicated.
            logger.warn("GVCF format is currently incompatible with allele validation. Not validating Alleles.");
            excludeTypesTemp.add(ValidationType.ALLELES);
        }
        final Set<ValidationType> excludeTypeSet = new LinkedHashSet<>(excludeTypesTemp);
        if (excludeTypesTemp.size() != excludeTypeSet.size()) {
            logger.warn("found repeat redundant validation types listed using the --validation-type-to-exclude argument");
        }
        if (excludeTypeSet.contains(ValidationType.ALL)) {
            if (excludeTypeSet.size() > 1) {
                logger.warn("found ALL in the --validation-type-to-exclude list together with other concrete type exclusions that are redundant");
            }
            final boolean[] allFalse = new boolean[ValidationType.CONCRETE_TYPES.size()];
            for (ValidationType t : ValidationType.CONCRETE_TYPES) {
                allFalse[t.ordinal()] = false;
            }
            return allFalse;
        } else {
            final Set<ValidationType> result = new LinkedHashSet<>(ValidationType.DEFAULT_TYPES);
            result.removeAll(excludeTypeSet);
            includeTypes = result;
        }
        includeTypes.addAll(this.includeTypes);
        if (this.includeTypes.contains(ValidationType.REF) && !hasReference()) {
            throw new UserException.MissingReference("Validation type " + ValidationType.REF.name() + " was selected but no reference was provided.");
        }
        if (includeTypes.contains(ValidationType.ALL)) {
            if (!hasReference()) {
                throw new UserException.MissingReference("Validation type " + ValidationType.REF.name() + " was selected but no reference was provided.");
            }
            final boolean[] allTrues = new boolean[ValidationType.CONCRETE_TYPES.size()];
            for (ValidationType t : ValidationType.CONCRETE_TYPES) {
                allTrues[t.ordinal()] = true;
            }
            return allTrues;
        }
        final boolean[] returnArray = new boolean[ValidationType.values().length];
        for (ValidationType t : includeTypes) {
            returnArray[t.ordinal()] = true;
        }
        return returnArray;
    }

    private void validateVariantsOrder(final VariantContext vc) {
        // Check if the current VC belongs to the same contig as the previous one.
        // If not, reset the start position to -1.
        if (previousContig == null || !previousContig.equals(vc.getContig())) {
            previousContig = vc.getContig();
            previousStart = -1;
        }

        //if next VC refers to a previous genomic position, throw an error
        //Note that HaplotypeCaller can emit variants that start inside of a deletion on another haplotype,
        // making v2's start less than the deletion's end
        if (previousStart > -1 && vc.getStart() < previousStart) {
            final UserException e = new UserException(String.format("In a GVCF all records must ordered. Record: %s covers a position previously traversed.",
                    vc.toStringWithoutGenotypes()));
            throwOrWarn(e);
        }
    }

    private void validateGVCFVariant(final VariantContext vc) {
        if (!vc.hasAllele(Allele.NON_REF_ALLELE)) {
            final UserException e = new UserException(String.format("In a GVCF all records must contain a %s allele. Offending record: %s",
                    Allele.NON_REF_STRING, vc.toStringWithoutGenotypes()));
            throwOrWarn(e);
        }
    }

    /**
     * May be run on GVCFs to ensure that single-sample pipeline outputs have required data for genotyping and filtering
     * @param vc
     */
    private void validateVQSRInputs(final VariantContext vc) {
        if (vc.isReferenceBlock() && !VALIDATE_GVCF) {
            final UserException e = new UserException.BadInput("VQSR should not be run on GVCFs.  Encountered a reference block at " +
                    vc.getContig() + ":" + vc.getStart() + " during VQSR input validation.");
            throwOrWarn(e);
        }
        if (vc.isFiltered()) {
            return;
        }
        if (!EXOME_INPUT) {
            checkForAnnotation(vc, VCFConstants.DEPTH_KEY);
        }
        boolean hasHetCall = false;
        if (vc.hasGenotypes()) {
            for (final Genotype g : vc.getGenotypes()) {
                if (needsRankSum(g)) {
                    hasHetCall = true;
                    break;
                }
            }
        } else {
            oneShotLogger.warn("No genotypes are present for variant at "
                    + vc.getContig() + ":" + vc.getStart() + "-- will not validate the existence of RankSum annotations");
        }
        if (VALIDATE_GVCF) {
            validateRequiredRawVQSRAnnotations(vc, hasHetCall);
            if (validationsToPerform[ValidationType.AS_ANNOTATIONS.ordinal()]) {
                validateRequiredRawASVQSRAnnotations(vc, hasHetCall);
            }
        } else {
            validateRequiredVQSRAnnotations(vc, hasHetCall);
            if (validationsToPerform[ValidationType.AS_ANNOTATIONS.ordinal()]) {
                validateRequiredASVQSRAnnotations(vc, hasHetCall);
            }
        }
    }

    private void applyValidations(final VariantContext vc, final Allele reportedRefAllele, final Allele observedRefAllele, final Set<String> rsIDs) {
        // Note: VariantContext.validateRSIDs blows up on an empty list (but works fine with null).
        // The workaround is to not pass an empty list.
       ValidationType currentType = null;

       try {
           if (validationsToPerform[ValidationType.REF.ordinal()]) {
               if (hasReference()) {
                   currentType = ValidationType.REF;
                   vc.validateReferenceBases(reportedRefAllele, observedRefAllele);
               }
           }
           if (validationsToPerform[ValidationType.IDS.ordinal()]) {
               if (!rsIDs.isEmpty()) {
                   currentType = ValidationType.IDS;
                   vc.validateRSIDs(rsIDs);
               }
           }
           if (validationsToPerform[ValidationType.ALLELES.ordinal()]) {
               currentType = ValidationType.ALLELES;
               vc.validateAlternateAlleles();
           }
           if (validationsToPerform[ValidationType.CHR_COUNTS.ordinal()]) {
               currentType = ValidationType.CHR_COUNTS;
               vc.validateChromosomeCounts();
           }
           if (validationsToPerform[ValidationType.GNARLY_INPUT.ordinal()]) {
               currentType = ValidationType.GNARLY_INPUT;
               validateGnarlyInputs(vc);
           }
           if (validationsToPerform[ValidationType.VQSR_INPUT.ordinal()]) {
               currentType = ValidationType.VQSR_INPUT;
               if (VALIDATE_GVCF) {
                   validateRequiredRawVQSRAnnotations(vc, needsRankSum(vc.getGenotype(0))  && vc.getGenotype(0).getAlleles().stream().anyMatch(Allele::isReference));
                   if (validationsToPerform[ValidationType.AS_ANNOTATIONS.ordinal()]) {
                       validateRequiredRawASVQSRAnnotations(vc, needsRankSum(vc.getGenotype(0)) && vc.getGenotype(0).getAlleles().stream().anyMatch(Allele::isReference));
                   }
               } else {
                   validateVQSRInputs(vc);
                   if (validationsToPerform[ValidationType.AS_ANNOTATIONS.ordinal()]) {
                       validateRequiredASVQSRAnnotations(vc, needsRankSum(vc.getGenotype(0))  && vc.getGenotype(0).getAlleles().stream().anyMatch(Allele::isReference));
                   }
               }
           }
           if (validationsToPerform[ValidationType.AS_ANNOTATIONS.ordinal()]) {
               currentType = ValidationType.AS_ANNOTATIONS;
               validateAlleleSpecificAnnotations(vc);
           }
       } catch (final UserException e) {
           final UserException withVC = new UserException.FailsStrictValidation(drivingVariantFile.toString(), currentType,
                   e.getMessage() + "Failure at variant context: " + vc.toStringWithoutGenotypes());
           throwOrWarn(withVC);
       }
    }

    private void validateGnarlyInputs(final VariantContext vc) {
        if (VALIDATE_GVCF && vc.getGenotype(0).isHomRef()) {
            return;
        }
        for (final String ann : GnarlyGenotyper.GNARLY_REQUIRED_INPUT_ANNOTATIONS) {
            checkForAnnotation(vc, ann);
        }
        if (validationsToPerform[ValidationType.AS_ANNOTATIONS.ordinal()]) {
            for (final String ann : GnarlyGenotyper.GNARLY_REQUIRED_AS_INPUT_ANNOTATIONS) {
                checkForAnnotation(vc, ann);
            }
        }
    }

    private void validateGnarlyAnnotations(final VariantContext vc) {
        for (final String ann : GnarlyGenotyper.GNARLY_EXPECTED_OUTPUT_ANNOTATIONS) {
            checkForAnnotation(vc, ann);
        }
        if (validationsToPerform[ValidationType.AS_ANNOTATIONS.ordinal()]) {
            for (final String ann : GnarlyGenotyper.GNARLY_EXPECTED_AS_OUTPUT_ANNOTATIONS) {
                checkForAnnotation(vc, ann);
            }
        }
    }

    private void checkForAnnotation(final VariantContext vc, final String annotationKey) {
        if (!vc.hasAttribute(annotationKey)) {
            final UserException e = new UserException.BadInput("Variant at " + vc.getContig() + ":" + vc.getStart() + " is missing " + annotationKey);
            throwOrWarn(e);
        }
    }

    private boolean needsRankSum(final Genotype g) {
        return g.isHet() && !g.isHetNonRef() && g.getAD()[0] > 0;
    }

    private void validateGnarlyOutputs(final VariantContext vc) {
        validateGnarlyAnnotations(vc);
        validateGnarlyGenotypes(vc);
    }

    private void validateGnarlyGenotypes(final VariantContext vc) {
        final int nObservedAlts = vc.getNAlleles() - 1;
        final boolean shouldHavePLs = nObservedAlts <= MAX_ALT_ALLELES;
        for (final Genotype g : vc.getGenotypes()) {
            if (g.isNoCall()) {
                break;
            }
            if (!g.isHomRef()) {
                validateRequiredGnarlyVariantGTAttributes(g, vc);
            }
            if (shouldHavePLs) {
                validateAdditionalGnarlyGTAttributes(g, vc);
            }
        }
    }

    /**
     * Check that variant genotypes have appropriate annotations
     * @param g
     * @param vc parent variant context, for error reporting
     */
    private void validateRequiredGnarlyVariantGTAttributes(final Genotype g, final VariantContext vc) {
        if (!g.hasAD()) {
            if (g.countAllele(Allele.SPAN_DEL) + g.countAllele(vc.getReference()) < 2) {
                final UserException e = new UserException.BadInput("Genotype for sample " + g.getSampleName() + " is missing AD at " + vc.getContig() + ":" + vc.getStart() + " : " + g);
                throwOrWarn(e);
            }
        }
    }

    private void validateAdditionalGnarlyGTAttributes(final Genotype g, final VariantContext vc) {
        if (!g.hasPL()) {
            final UserException e = new UserException.BadInput("Genotype for sample " + g.getSampleName() + " is missing PL at " + vc.getContig() + ":" + vc.getStart() + " : " + g);
            throwOrWarn(e);
        }
        if (!g.hasGQ()) {
            final UserException e = new UserException.BadInput("Genotype for sample " + g.getSampleName() + " is missing GQ at " + vc.getContig() + ":" + vc.getStart() + " : " + g);
            throwOrWarn(e);
        }
    }

    private void validateAlleleSpecificAnnotations(final VariantContext vc) {
        final UserException e = GATKVariantContextUtils.assertAlleleSpecificAnnotationsHaveCorrectLength(vc);
        if (e != null) {
            throwOrWarn(e);
        }
    }

    private void validateRequiredVQSRAnnotations(final VariantContext vc, final boolean hasHetCall) {
        for (final String requiredAnnotation : requiredVQSRAnnotationKeys) {
            if (requiredAnnotation.contains("RankSum") && !hasHetCall) {
                continue;
            }
            checkForAnnotation(vc, requiredAnnotation);
        }
    }

    private void validateRequiredRawVQSRAnnotations(final VariantContext vc, final boolean hasHetCall) {
        if (VALIDATE_GVCF && vc.getGenotype(0).isHomRef()) {
            return;
        }
        if (!vc.hasGenotypes()) {
            throw new UserException("Strand bias annotations cannot be calculated without FORMAT-level SB field.");
        }
        for (final Genotype gt : vc.getGenotypes()) {
            if (!gt.hasAnyAttribute(GATKVCFConstants.STRAND_BIAS_BY_SAMPLE_KEY)) {
                final UserException e = new UserException.BadInput("Sample " + gt.getSampleName() + " in variant at " +
                        vc.getContig() + ":" + vc.getStart() + " is missing " + GATKVCFConstants.STRAND_BIAS_BY_SAMPLE_KEY);
                throwOrWarn(e);
            }
        }
        for (final String requiredAnnotation : requiredRawVQSRAnnotationKeys) {
            if (requiredAnnotation.contains("RankSum") && !hasHetCall) {
                continue;
            }
            checkForAnnotation(vc, requiredAnnotation);
        }
    }

    private void validateRequiredASVQSRAnnotations(final VariantContext vc, final boolean hasHetCall) {
        for (final String requiredAnnotation : requiredAlleleSpecificVQSRAnnotationKeys) {
            if (requiredAnnotation.contains("RankSum") && !hasHetCall) {
                continue;
            }
            checkForAnnotation(vc, requiredAnnotation);
        }
    }

    private void validateRequiredRawASVQSRAnnotations(final VariantContext vc, final boolean hasHetCall) {
        if (VALIDATE_GVCF && vc.getGenotype(0).isHomRef()) {
            return;
        }
        for (final String requiredAnnotation : requiredRawASVQSRAnnotationKeys) {
            if (requiredAnnotation.contains("RankSum") && !hasHetCall) {
                continue;
            }
            checkForAnnotation(vc, requiredAnnotation);
        }
    }

    private void throwOrWarn(UserException e) {
        if (WARN_ON_ERROR) {
            logger.warn("***** " + e.getMessage() + " *****");
        } else {
            throw e;
        }
    }
}
