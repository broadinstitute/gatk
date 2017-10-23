package org.broadinstitute.hellbender.tools.walkers.validation.basicshortmutpileup;


import htsjdk.variant.variantcontext.Allele;
import org.apache.commons.math3.distribution.BinomialDistribution;
import org.broadinstitute.hellbender.utils.GATKProtectedVariantContextUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.pileup.PileupElement;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;

import java.util.List;
import java.util.stream.Collectors;

public class PowerCalculationUtils {
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

    static int calculateMinCountForSignal(int validationTumorTotalCount, double maxSignalRatioInNormal) {
        final BinomialDistribution binomialDistribution = new BinomialDistribution(validationTumorTotalCount, maxSignalRatioInNormal);
        int minCountForSignal = binomialDistribution.inverseCumulativeProbability(0.99);
        minCountForSignal = Math.max(minCountForSignal, 2);
        return minCountForSignal;
    }

    /**
     * Calculate the maximum ratio of (non-ref read count) to (reference allele read count) + (non-ref read count)
     *  in the given pileup.
     *
     *  In this context, this code is only looking for reads that have a base change (or directly precede an indel) at that location.
     *
     * @param readPileup pileup to create the ratio.
     * @param referenceAllele reference allele
     * @param minBaseQualityCutoff only count the bases that exceed the min base quality.  For xNP, it is the min quality
     *                             all overlapping loci.  For indels, this is the base preceding the indel itself.
     * @return ratio as a double
     */
    static double calculateMaxAltRatio(final ReadPileup readPileup, final Allele referenceAllele, int minBaseQualityCutoff) {

        // Go through the read pileup and find all new bases.
        final List<PileupElement> pileupElementsPassingQuality = Utils.stream(readPileup.iterator())
                .filter(pe -> pe.getBase() != PileupElement.DELETION_BASE)
                .filter(pe -> pe.getRead().getBaseQuality(pe.getOffset()) >= minBaseQualityCutoff)
                .collect(Collectors.toList());

        final long numAlternate = pileupElementsPassingQuality.stream()
                .filter(pe -> (GATKProtectedVariantContextUtils.doesReadContainAllele(pe, referenceAllele) == GATKProtectedVariantContextUtils.ReadContainAllele.FALSE)
                || pe.isBeforeDeletionStart() || pe.isBeforeInsertion()).count();

        final long numReference = pileupElementsPassingQuality.stream()
                .filter(pe -> (GATKProtectedVariantContextUtils.doesReadContainAllele(pe, referenceAllele) == GATKProtectedVariantContextUtils.ReadContainAllele.TRUE)
                && !pe.isBeforeDeletionStart() && !pe.isBeforeInsertion()).count();

        return (double) numAlternate / ((double) numReference + (double) numAlternate);
    }
}
