package org.broadinstitute.hellbender.tools.walkers.variantutils;

import htsjdk.tribble.TribbleException;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.ArgumentCollection;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.argumentcollections.DbsnpArgumentCollection;
import org.broadinstitute.hellbender.cmdline.programgroups.VariantProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.VariantWalker;
import org.broadinstitute.hellbender.exceptions.UserException;

import java.util.*;


/**
 * Validates a VCF file with an extra strict set of criteria.
 *
 * <p>
 * ValidateVariants is a GATK tool that takes a VCF file and validates much of the information inside it.
 * In addition to standard adherence to the VCF specification, this tool performs extra strict validations to ensure
 * the information contained within the file is correct. These include:
 * </p><p>
 * <dl>
 *   <dt>REF</dt><dd>the correctness of the reference base(s).</dd>
 *   <dt>CHR_COUNTS</dt><dd>accuracy of AC & AN values.</dd>
 *   <dt>IDS</dt><dd>tests against rsIDs when a dbSNP file is provided. Notice that for this one to work, you need
 *    to provide a reference to the dbsnp variant containing file using the <code>--dbsnp</code> as show in examples below.</dd>
 *   <dt>ALLELES</dt><dd>and that all alternate alleles are present in at least one sample.</dd>
 * </dl>
 *
 * </p>
 *
 * <p>
 *     By default it will apply all the strict validations unless you indicate which one you want you want to exclude
 *     using <code>-Xtype|--validationTypeToExclude &lt;<i>code</i>&gt;</code>, where <i>code</i> is one of the listed above. You
 *     can exclude as many types as you want
 * <p>
 *     You can exclude all strict validations with the special code <code><b>ALL</b></code>. In this case the tool will only
 *     test the adherence to the VCF specification.
 * </p>
 *
 * <h3>Input</h3>
 * <p>
 * A variant set to validate using <code>-V</code> or <code>--variant</code> as shown below.
 * </p>
 *
 * <h3>Examples</h3>
 *
 * <p>To perform VCF format and all strict validations: </p>
 *
 * <pre>
 * java -Xmx2g -jar GenomeAnalysisTK.jar \
 *   -R ref.fasta \
 *   -T ValidateVariants \
 *   --variant input.vcf \
 *   --dbsnp dbsnp.vcf
 * </pre>
 *
 * <p>To perform only VCF format tests:</p>
 *
 * <pre>
 * java -Xmx2g -jar GenomeAnalysisTK.jar \
 *   -R ref.fasta \
 *   -T ValidateVariants \
 *   <b>--validationTypeToExclude ALL</b> \
 *   --variant input.vcf
 * </pre>
 *
 * <p>To perform all validations except the strict <i>ALLELE</i> validation:</p>
 *
 * <pre>
 * java -Xmx2g -jar GenomeAnalysisTK.jar \
 *   -R ref.fasta \
 *   -T ValidateVariants \
 *   <b>--validationTypeToExclude ALLELES</b>
 *   --variant input.vcf \
 *   --dbsnp dbsnp.vcf
 * </pre>
 *
 */
@CommandLineProgramProperties(
        summary = "Validates a VCF file with an extra strict set of criteria.",
        oneLineSummary = "Validate VCF",
        programGroup = VariantProgramGroup.class
)
public final class ValidateVariants extends VariantWalker {
    static final Logger logger = LogManager.getLogger(ValidateVariants.class);

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

    @Argument(fullName = "validationTypeToExclude", shortName = "Xtype", doc = "which validation type to exclude from a full strict validation", optional = true)
    List<ValidationType> excludeTypes = new ArrayList<>();

    /**
     * By default, even filtered records are validated.
     */
    @Argument(fullName = "doNotValidateFilteredRecords", shortName = "doNotValidateFilteredRecords", doc = "skip validation on filtered records", optional = true)
    Boolean DO_NOT_VALIDATE_FILTERED = false;

    @Argument(fullName = "warnOnErrors", shortName = "warnOnErrors", doc = "just emit warnings on errors instead of terminating the run at the first instance", optional = true)
    Boolean WARN_ON_ERROR = false;

    /**
     * Contains final set of validation to apply.
     */
    private Collection<ValidationType> validationTypes;

    @Override
    public void onTraversalStart() {
        validationTypes = calculateValidationTypesToApply(excludeTypes);
    }

    @Override
    public void apply(final VariantContext vc, final ReadsContext readsContext, final ReferenceContext ref, final FeatureContext featureContext) {
        if ( DO_NOT_VALIDATE_FILTERED && vc.isFiltered() ) {
            return;
        }

        // get the true reference allele
        final Allele reportedRefAllele = vc.getReference();
        final int refLength = reportedRefAllele.length();

        final byte[] observedRefBases = Arrays.copyOf(ref.getBases(), refLength);
        final Allele observedRefAllele = Allele.create(observedRefBases);
        final Set<String> rsIDs = getRSIDs(featureContext);

        for (final ValidationType t : validationTypes) {
            try{
                applyValidationType(vc, reportedRefAllele, observedRefAllele, rsIDs, t);
            } catch (TribbleException e) {
                if ( WARN_ON_ERROR ) {
                    logger.warn("***** " + e.getMessage() + " *****");
                } else {
                    throw new UserException.FailsStrictValidation(drivingVariantFile, t, e.getMessage());
                }
            }
        }
    }

    /*
     *  Returns the list of RSIDs overlapping the current variant that we're walking over.
     *  If there's no RSID or if there was not dbsnp file passed in as an argument,
     *  an empty set is returned (and then no validation is performed, see applyValidationType.
     */
    private Set<String> getRSIDs(FeatureContext featureContext) {
        Set<String> rsIDs = new HashSet<>();
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
        if (excludeTypes.size() == 0) {
            return Collections.singleton(ValidationType.ALL);
        }
        final Set<ValidationType> excludeTypeSet = new LinkedHashSet<>(excludeTypes);
        if (excludeTypes.size() != excludeTypeSet.size()) {
            logger.warn("found repeat redundant validation types listed using the --validationTypeToExclude argument");
        }
        if (excludeTypeSet.contains(ValidationType.ALL)) {
            if (excludeTypeSet.size() > 1) {
                logger.warn("found ALL in the --validationTypeToExclude list together with other concrete type exclusions that are redundant");
            }
            return Collections.emptyList();
        } else {
           final Set<ValidationType> result = new LinkedHashSet<>(ValidationType.CONCRETE_TYPES);
           result.removeAll(excludeTypeSet);
           return result;
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
}