package org.broadinstitute.hellbender.tools.walkers.variantutils;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.tribble.TribbleException;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.argumentcollections.DbsnpArgumentCollection;
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

    public static final String GVCF_VALIDATE = "validate-GVCF";
    public static final String DO_NOT_VALIDATE_FILTERED_RECORDS = "do-not-validate-filtered-records";

    public enum ValidationType {

        /**
         * Makes reference to all extra-strict tests listed below.
         */
        ALL,

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
        CHR_COUNTS;

        /**
         * Unmodifiable set of concrete validation types.
         *
         * <p>These are all types except {@link #ALL}.</p>
         */
        public static final Set<ValidationType> CONCRETE_TYPES;

        static {
            final Set<ValidationType> cts = new LinkedHashSet<>(values().length - 1);
            for (final ValidationType v : values()) {
                if (v != ALL)
                    cts.add(v);
            }
            CONCRETE_TYPES = Collections.unmodifiableSet(cts);
        }
    }

    @ArgumentCollection
    DbsnpArgumentCollection dbsnp = new DbsnpArgumentCollection();

    @Argument(fullName = "validation-type-to-exclude",
            shortName = "Xtype",
            doc = "which validation type to exclude from a full strict validation",
            optional = true)
    List<ValidationType> excludeTypes = new ArrayList<>();

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
     * Contains final set of validation to apply.
     */
    private Collection<ValidationType> validationTypes;

    private GenomeLocSortedSet genomeLocSortedSet;

    // information to keep track of when validating a GVCF
    private SimpleInterval previousInterval;

    @Override
    public void onTraversalStart() {
        if (VALIDATE_GVCF) {
            final SAMSequenceDictionary seqDictionary = getBestAvailableSequenceDictionary();

            if (seqDictionary == null)
                throw new UserException("Validating a GVCF requires a sequence dictionary but no dictionary was able to be constructed from your input.");

            genomeLocSortedSet = new GenomeLocSortedSet(new GenomeLocParser(seqDictionary));
        }
        validationTypes = calculateValidationTypesToApply(excludeTypes);
    }

    @Override
    public void apply(final VariantContext vc, final ReadsContext readsContext, final ReferenceContext ref, final FeatureContext featureContext) {
        if (DO_NOT_VALIDATE_FILTERED && vc.isFiltered()) {
            return;
        }
        // get the true reference allele
        final Allele reportedRefAllele = vc.getReference();
        final int refLength = reportedRefAllele.length();

        final Allele observedRefAllele = hasReference() ? Allele.create(Arrays.copyOf(ref.getBases(), refLength)) : null;

        final Set<String> rsIDs = getRSIDs(featureContext);

        if (VALIDATE_GVCF) {
            final SimpleInterval refInterval = ref.getInterval();

            // GenomeLocSortedSet will automatically merge intervals that are overlapping when setting `mergeIfIntervalOverlaps`
            // to true.  In a GVCF most blocks are adjacent to each other so they wouldn't normally get merged.  We check
            // if the current record is adjacent to the previous record and "overlap" them if they are so our set is as
            // small as possible while still containing the same bases.
            final int start = (previousInterval != null && previousInterval.overlapsWithMargin(refInterval, 1)) ?
                    previousInterval.getStart() : refInterval.getStart();
            final int end = (previousInterval != null && previousInterval.overlapsWithMargin(refInterval, 1)) ?
                    Math.max(previousInterval.getEnd(), vc.getEnd()) : vc.getEnd();
            final GenomeLoc possiblyMergedGenomeLoc = genomeLocSortedSet.getGenomeLocParser().createGenomeLoc(refInterval.getContig(), start, end);
            genomeLocSortedSet.add(possiblyMergedGenomeLoc, true);

            previousInterval = new SimpleInterval(possiblyMergedGenomeLoc);
            validateGVCFVariant(vc);
        }

        for (final ValidationType t : validationTypes) {
            try{
                applyValidationType(vc, reportedRefAllele, observedRefAllele, rsIDs, t);
            } catch (TribbleException e) {
                throwOrWarn(new UserException.FailsStrictValidation(drivingVariantFile, t, e.getMessage()));
            }
        }
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
                final UserException e = new UserException("A GVCF must cover the entire region. Found " + uncoveredIntervals.coveredSize() +
                        " loci with no VariantContext covering it. The first uncovered segment is:" +
                        uncoveredIntervals.iterator().next());
                throwOrWarn(e);
            }
        }
        return null;
    }

    /*
     *  Returns the list of RSIDs overlapping the current variant that we're walking over.
     *  If there's no RSID or if there was not dbsnp file passed in as an argument,
     *  an empty set is returned (and then no validation is performed, see applyValidationType.
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
    private Collection<ValidationType> calculateValidationTypesToApply(final List<ValidationType> excludeTypes) {
        if (VALIDATE_GVCF && !excludeTypes.contains(ValidationType.ALLELES)) {
            // Note: in a future version allele validation might be OK for GVCFs, if that happens
            // this will be more complicated.
            logger.warn("GVCF format is currently incompatible with allele validation. Not validating Alleles.");
            excludeTypes.add(ValidationType.ALLELES);
        }
        if (excludeTypes.isEmpty()) {
            return Collections.singleton(ValidationType.ALL);
        }
        final Set<ValidationType> excludeTypeSet = new LinkedHashSet<>(excludeTypes);
        if (excludeTypes.size() != excludeTypeSet.size()) {
            logger.warn("found repeat redundant validation types listed using the --validation-type-to-exclude argument");
        }
        if (excludeTypeSet.contains(ValidationType.ALL)) {
            if (excludeTypeSet.size() > 1) {
                logger.warn("found ALL in the --validation-type-to-exclude list together with other concrete type exclusions that are redundant");
            }
            return Collections.emptyList();
        } else {
            final Set<ValidationType> result = new LinkedHashSet<>(ValidationType.CONCRETE_TYPES);
            result.removeAll(excludeTypeSet);
            if (result.contains(ValidationType.REF) && !hasReference()) {
                throw new UserException.MissingReference("Validation type " + ValidationType.REF.name() + " was selected but no reference was provided.");
            }
            return result;
        }
    }

    private void validateGVCFVariant(final VariantContext vc) {
        if (!vc.hasAllele(Allele.NON_REF_ALLELE)) {
            final UserException e = new UserException(String.format("In a GVCF all records must contain a %s allele. Offending record: %s",
                    GATKVCFConstants.NON_REF_SYMBOLIC_ALLELE_NAME, vc.toStringWithoutGenotypes()));
            throwOrWarn(e);
        }
    }

    private void applyValidationType(VariantContext vc, Allele reportedRefAllele, Allele observedRefAllele, Set<String> rsIDs, ValidationType t) {
        // Note: VariantContext.validateRSIDs blows up on an empty list (but works fine with null).
        // The workaround is to not pass an empty list.
        switch( t ) {
            case ALL:
                if (!rsIDs.isEmpty()) {
                    vc.extraStrictValidation(reportedRefAllele, observedRefAllele, rsIDs);
                }
                break;
            case REF:
                vc.validateReferenceBases(reportedRefAllele, observedRefAllele);
                break;
            case IDS:
                if (!rsIDs.isEmpty()) {
                    vc.validateRSIDs(rsIDs);
                }
                break;
            case ALLELES:
                vc.validateAlternateAlleles();
                break;
            case CHR_COUNTS:
                vc.validateChromosomeCounts();
                break;
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
