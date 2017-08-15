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
import org.broadinstitute.hellbender.cmdline.programgroups.VariantProgramGroup;
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
 * <p>This tool is designed to validate the correctness of the formatting of VCF files. In addition to standard
 * adherence to the VCF specification, this tool performs extra strict validations to ensure that the information
 * contained within the file is correctly encoded. These include:
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
 *     By default the tool applies all the strict validations unless you indicate which one you want you want to
 *     exclude using `--validationTypeToExclude`. You can exclude as many types as you want. You can exclude all strict
 *     validations with the special code `ALL`. In this case the tool will only test for adherence to the VCF
 *     specification.
 * </p>
 *
 * <h3>Input</h3>
 * <p>
 * A VCF file to validate.
 * </p>
 *
 * <h3>Usage examples</h3>
 *
 * <h4>To perform VCF format and all strict validations: </h4>
 * <pre>
 *   ./gatk-launch ValidateVariants \
 *   -R ref.fasta \
 *   -V input.vcf \
 *   --dbsnp dbsnp.vcf
 * </pre>
 *
 * <h4>To perform only VCF format tests:</h4>
 * <pre>
 *   ./gatk-launch ValidateVariants
 *   -R ref.fasta \
 *   -V input.vcf \
 *   --validationTypeToExclude ALL
 * </pre>
 *
 * <h4>To perform all validations except the strict `ALLELE` validation:</h4>
 * <pre>
 *   ./gatk-launch ValidateVariants \
 *   -R ref.fasta \
 *   -V input.vcf \
 *   --validationTypeToExclude ALLELES \
 *   --dbsnp dbsnp.vcf
 * </pre>
 *
 */
@CommandLineProgramProperties(
        summary = "Validates a VCF file with an extra strict set of criteria.",
        oneLineSummary = "Validate VCF",
        programGroup = VariantProgramGroup.class
)
@DocumentedFeature
public final class ValidateVariants extends VariantWalker {
    static final Logger logger = LogManager.getLogger(ValidateVariants.class);

    public static final String GVCF_VALIDATE = "validateGVCF";
    public static final String DO_NOT_VALIDATE_FILTERED_RECORDS = "doNotValidateFilteredRecords";

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
    @Argument(fullName = DO_NOT_VALIDATE_FILTERED_RECORDS, shortName = "doNotValidateFilteredRecords", doc = "skip validation on filtered records", optional = true, mutex = GVCF_VALIDATE)
    Boolean DO_NOT_VALIDATE_FILTERED = false;

    @Argument(fullName = "warnOnErrors", shortName = "warnOnErrors", doc = "just emit warnings on errors instead of terminating the run at the first instance", optional = true)
    Boolean WARN_ON_ERROR = false;

    /**
     *  Validate this file as a gvcf. In particular, every variant must include a <NON_REF> allele, and that
     *  every base in the territory under consideration is covered by a variant (or a reference block).
     *  If you specifed intervals (using -L or -XL) to restrict analysis to a subset of genomic regions,
     *  those intervals will need to be covered in a valid gvcf.
     */
    @Argument(fullName = GVCF_VALIDATE, shortName = "gvcf", doc = "Validate this file as a GVCF", optional = true, mutex = DO_NOT_VALIDATE_FILTERED_RECORDS)
    Boolean VALIDATE_GVCF = false;

    /**
     * Contains final set of validation to apply.
     */
    private Collection<ValidationType> validationTypes;

    private GenomeLocSortedSet genomeLocSortedSet;

    // information to keep track of when validating a GVCF
    private int currentEnd;
    private String currentContig;

    SAMSequenceDictionary dict;

    @Override
    public void onTraversalStart() {
        if (VALIDATE_GVCF) {
            dict = getBestAvailableSequenceDictionary();

            if (dict == null)
                throw new UserException("Validating a GVCF requires a sequence dictionary but no dictionary was able to be constructed from your input.");

            genomeLocSortedSet = new GenomeLocSortedSet(new GenomeLocParser(dict));
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
            SimpleInterval refInterval = ref.getInterval();
            // Take advantage of adjacent blocks and just merge them so we dont have to keep so many objects in the set.
            if (refInterval.getContig().equals(currentContig) && currentEnd == (refInterval.getStart() - 1)) {
                genomeLocSortedSet.add(genomeLocSortedSet.getGenomeLocParser().createGenomeLoc(refInterval.getContig(), currentEnd, vc.getEnd()), true);
            } else {
                genomeLocSortedSet.add(genomeLocSortedSet.getGenomeLocParser().createGenomeLoc(refInterval.getContig(), refInterval.getStart(), vc.getEnd()), true);
            }
            currentEnd = vc.getEnd();
            currentContig = vc.getContig();
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
            GenomeLocSortedSet intervalArgumentGenomeLocSortedSet;

            if (intervalArgumentCollection.intervalsSpecified()){
                intervalArgumentGenomeLocSortedSet = GenomeLocSortedSet.createSetFromList(genomeLocSortedSet.getGenomeLocParser(), IntervalUtils.genomeLocsFromLocatables(genomeLocSortedSet.getGenomeLocParser(), intervalArgumentCollection.getIntervals(dict)));
            } else {
                intervalArgumentGenomeLocSortedSet = GenomeLocSortedSet.createSetFromSequenceDictionary(dict);
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
            if (result.contains(ValidationType.REF) && !hasReference()) {
                throw new UserException.MissingReference("Validation type " + ValidationType.REF.name() + " was selected but no reference was provided.");
            }
           return result;
        }
    }

    private void validateGVCFVariant(final VariantContext vc) {
        if (!vc.hasAllele(GATKVCFConstants.NON_REF_SYMBOLIC_ALLELE)) {
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
