package org.broadinstitute.hellbender.tools.walkers.validation.basicshortmutpileup;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hdf5.Utils;
import org.broadinstitute.hellbender.utils.GATKProtectedVariantContextUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;

public class BasicSomaticShortMutationValidator {
    private static final Logger logger = LogManager.getLogger(BasicSomaticShortMutationValidator.class);

    /**
     * @param genotype The genotype to test whether we can attempt to validate.  Never {@code null}
     * @param referenceAllele The reference allele (from parent VariantContext).  Never {@code null}
     * @return whether this class can even attempt a validation of the genotype in question.
     */
    public static boolean isAbleToValidateGenotype(final Genotype genotype, final Allele referenceAllele) {
        Utils.nonNull(genotype);
        Utils.nonNull(referenceAllele);

        // In order to proceed, we have some assumptions:
        //  - The genotype has a ploidy of 2
        //  - The genotype is a simple indel or a xNP
        //  - The first allele of the genotype is reference
        final boolean isPloidyOfTwo = genotype.getAlleles().size() == 2;
        final boolean doesGenotypeHaveReference = genotype.getAllele(0).equals(referenceAllele);
        final VariantContext.Type typeOfVariant = GATKProtectedVariantContextUtils.typeOfVariant(genotype.getAllele(0), genotype.getAllele(1));
        final boolean isValidatableVariantType =
                (typeOfVariant.equals(VariantContext.Type.INDEL) || typeOfVariant.equals(VariantContext.Type.SNP) ||
                        typeOfVariant.equals(VariantContext.Type.MNP))
                        && !GATKProtectedVariantContextUtils.isComplexIndel(genotype.getAllele(0), genotype.getAllele(1));
        final boolean hasKnownCoverage = genotype.hasAD() && (genotype.getAD().length == 2);
        final boolean isCanValidate = (isPloidyOfTwo && doesGenotypeHaveReference && isValidatableVariantType && hasKnownCoverage);
        if (!isCanValidate) {
            logger.info("Cannot validate genotype: " + genotype + "  ploidy2: " + isPloidyOfTwo +
                    "  genotypeHasReferenceAllele: " + doesGenotypeHaveReference + "   validatableVariant: " + isValidatableVariantType +
                    "  hasCompleteAD field: " + hasKnownCoverage);
        }

        return isCanValidate;
    }

    private BasicSomaticShortMutationValidator() {
    }

    /** Perform basic somatic pileup validation and return a result instance.
     *
     * @param genotype given genotype to validate.  Never {@code null}
     * @param referenceAllele Never {@code null}
     * @param validationNormalPileup Pileup for the coresponding location of the genotype in the validation normal.  Never {@code null}
     * @param validationTumorAltCount Alt read count for the validation tumor sample.  Must be positive or zero.
     * @param validationTumorTotalCount Total read count for the validation tumor sample.  Must be positive or zero.
     * @param minBaseQualityCutoff Minimum base quality threshold to count any read.  Must be positive or zero.
     * @param interval location represented in the genotype.  Never {@code null}
     * @param filters additional filters besides the ones on the genotype.  Typically,
     *                the parent variant context filters is provided here.  Never {@code null}
     * @return validation result.  {@code null} if the genotype cannot be validated
     */
    public static BasicValidationResult calculateBasicValidationResult(final Genotype genotype, final Allele referenceAllele,
                                                                       final ReadPileup validationNormalPileup,
                                                                       int validationTumorAltCount, int validationTumorTotalCount,
                                                                       int minBaseQualityCutoff, final SimpleInterval interval, final String filters) {

        if (!isAbleToValidateGenotype(genotype, referenceAllele)){
            return null;
        }

        Utils.nonNull(referenceAllele);
        Utils.nonNull(genotype);
        Utils.nonNull(validationNormalPileup);
        Utils.nonNull(interval);
        Utils.nonNull(filters);
        ParamUtils.isPositiveOrZero(validationTumorAltCount, "Validation alt count must be >= 0");
        ParamUtils.isPositiveOrZero(validationTumorTotalCount, "Validation total count must be >= 0");
        ParamUtils.isPositiveOrZero(minBaseQualityCutoff, "Minimum base quality cutoff must be >= 0");

        final double maxAltRatioSeenInNormalValidation = PowerCalculationUtils.calculateMaxAltRatio(validationNormalPileup,
                referenceAllele, minBaseQualityCutoff);
        final int discoveryTumorAltCount = genotype.getAD()[1];
        final int discoveryTumorTotalCount = genotype.getAD()[0] + discoveryTumorAltCount;

        final int minCountForSignal = PowerCalculationUtils.calculateMinCountForSignal(validationTumorTotalCount, maxAltRatioSeenInNormalValidation);
        final double power = PowerCalculationUtils.calculatePower(validationTumorTotalCount, discoveryTumorAltCount, discoveryTumorTotalCount, minCountForSignal);

        boolean isNotNoise = (validationTumorAltCount >= minCountForSignal);
        boolean isEnoughValidationCoverageToValidate = validationTumorAltCount >= 2;

        final String genotypeFilters = genotype.getFilters() == null ? "" : genotype.getFilters();

        return new BasicValidationResult(interval, minCountForSignal, isEnoughValidationCoverageToValidate,
                isNotNoise, power, validationTumorAltCount, validationTumorTotalCount-validationTumorAltCount,
                discoveryTumorAltCount, discoveryTumorTotalCount-discoveryTumorAltCount, referenceAllele,
                genotype.getAllele(1), filters + genotypeFilters);
    }
}
