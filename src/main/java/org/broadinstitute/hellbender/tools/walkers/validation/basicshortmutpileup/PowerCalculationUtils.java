package org.broadinstitute.hellbender.tools.walkers.validation.basicshortmutpileup;


import htsjdk.variant.variantcontext.Allele;
import org.apache.commons.math3.distribution.BinomialDistribution;
import org.broadinstitute.hellbender.utils.GATKProtectedVariantContextUtils;
import org.broadinstitute.hellbender.utils.Trilean;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.broadinstitute.hellbender.utils.pileup.PileupElement;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;

import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

public class PowerCalculationUtils {

    public static final double P_VALUE_FOR_NOISE = 0.99;
    public static final int MINIMUM_NUM_READS_FOR_SIGNAL_COUNT = 2;

    private PowerCalculationUtils() {}

    /**
     *  Calculate the power to validate the mutation with given alt and total counts.
     *
     * @param validationTumorTotalCount total read count in the validation tumor bam
     * @param discoveryTumorAltCount alt read count in the discovery tumor
     * @param discoveryTumorTotalCount total read count in the discovery tumor
     * @param minCountForSignal minimum number of alt reads to claim that the variant is out of the noise floor.
     * @return p-value that the discovery variant is out of the noise floor.
     */
    public static double calculatePower(int validationTumorTotalCount, int discoveryTumorAltCount, int discoveryTumorTotalCount, int minCountForSignal) {
        // n is the total count in the validation.  alpha is the alt count in the discovery tumor +1
        // beta is the ref count in the discovery tumor + 1
        final BetaBinomialDistribution betaBinomialDistribution =
                new BetaBinomialDistribution(null, discoveryTumorAltCount + 1, discoveryTumorTotalCount - discoveryTumorAltCount + 1, validationTumorTotalCount);
        return (1 - betaBinomialDistribution.cumulativeProbability(minCountForSignal - 1));
    }

    /**
     *
     * @param validationTumorTotalCount the total number of reads that we would like to validate.
     *                                  Must be positive or zero.
     * @param maxSignalRatioInNormal a ratio of non-ref to ref bases.  Must be within 0.0 - 1.0.
     * @return the minimum read count to be above the model noise floor.
     */
    public static int calculateMinCountForSignal(int validationTumorTotalCount, double maxSignalRatioInNormal) {
        ParamUtils.isPositiveOrZero(validationTumorTotalCount, "Cannot have a negative total count.");
        ParamUtils.inRange(maxSignalRatioInNormal, 0.0, 1.0, "Cannot have have a ratio that is outside of 0.0 - 1.0.");

        final BinomialDistribution binomialDistribution = new BinomialDistribution(validationTumorTotalCount, maxSignalRatioInNormal);
        return Math.max(binomialDistribution.inverseCumulativeProbability(P_VALUE_FOR_NOISE), MINIMUM_NUM_READS_FOR_SIGNAL_COUNT);
    }

    /**
     * Calculate the fraction of non-ref reads in the given pileup.
     *
     *  In this context, this code is only looking for reads that have a base change (or directly precede an indel) at that location.
     *
     * @param readPileup pileup to create the ratio.  Never {@code null}
     * @param referenceAllele reference allele.  Never {@code null}
     * @param minBaseQualityCutoff only count the bases that exceed the min base quality.  For xNP, it is the min quality
     *                             all overlapping loci.  For indels, this is the base preceding the indel itself.  Must be positive or zero.
     *                             Zero indicates that all reads should pass this filter.
     * @return ratio as a double
     */
    public static double calculateMaxAltRatio(final ReadPileup readPileup, final Allele referenceAllele, int minBaseQualityCutoff) {
        ParamUtils.isPositiveOrZero(minBaseQualityCutoff, "Cannot have a negative minBaseQualityCutoff.");
        Utils.nonNull(readPileup);
        Utils.nonNull(referenceAllele);

        // Go through the read pileup and find all new bases.
        final List<PileupElement> pileupElementsPassingQuality = retrievePileupElements(readPileup, minBaseQualityCutoff);

        final long numAlternate = pileupElementsPassingQuality.stream()
                .filter(pe -> (GATKProtectedVariantContextUtils.doesReadContainAllele(pe, referenceAllele) == Trilean.FALSE)
                || pe.isBeforeDeletionStart() || pe.isBeforeInsertion()).count();

        final long numReference = pileupElementsPassingQuality.stream()
                .filter(pe -> (GATKProtectedVariantContextUtils.doesReadContainAllele(pe, referenceAllele) == Trilean.TRUE)
                && !pe.isBeforeDeletionStart() && !pe.isBeforeInsertion()).count();

        return numReference + numAlternate == 0 ? 0.0 : (double) numAlternate / ((double) numReference + (double) numAlternate);
    }

    private static List<PileupElement> retrievePileupElements(final ReadPileup readPileup, final int minBaseQualityCutoff) {
        return Utils.stream(readPileup.iterator())
                    .filter(pe -> !pe.isDeletion())
                    .filter(pe -> pe.getQual() >= minBaseQualityCutoff)
                    .collect(Collectors.toList());
    }

    /**
     * Get a count of the number of reads supporting an allele in a given pileup.
     *
     * @param readPileup pileup to query for the allele.  Never {@code null}
     * @param referenceAllele  reference allele corresponding to alt allele.  Never {@code null}
     * @param altAllele Allele to query. Never {@code null}
     * @param minBaseQualityCutoff only count the bases that exceed the min base quality.  For xNP, it is the min quality
     *                             all overlapping loci.  For indels, this is the base preceding the indel itself.  Must be positive or zero.
     *                             Zero indicates that all reads should pass this filter.
     * @return count of reads supporting the given allele in the given pileup.
     */
    public static long calculateNumReadsSupportingAllele(final ReadPileup readPileup, final Allele referenceAllele, final Allele altAllele, int minBaseQualityCutoff) {
        ParamUtils.isPositiveOrZero(minBaseQualityCutoff, "Cannot have a negative minBaseQualityCutoff.");
        Utils.nonNull(readPileup);
        Utils.nonNull(altAllele);
        Utils.nonNull(referenceAllele);

        final List<PileupElement> pileupElementsPassingQuality = retrievePileupElements(readPileup, minBaseQualityCutoff);

        return pileupElementsPassingQuality.stream()
                .filter(pe -> altAllele.equals(GATKProtectedVariantContextUtils.chooseAlleleForRead(pe, referenceAllele, Collections.singletonList(altAllele), minBaseQualityCutoff)))
                .count();
    }
}
