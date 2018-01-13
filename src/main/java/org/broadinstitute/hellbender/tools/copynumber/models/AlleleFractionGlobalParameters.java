package org.broadinstitute.hellbender.tools.copynumber.models;

import org.broadinstitute.hellbender.tools.copynumber.formats.CopyNumberFormatsUtils;

/**
 * Encapsulates the global parameters of the allele fraction model: the mean and variance of the common prior on
 * allelic biases and the outlier probability.
 *
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 */
final class AlleleFractionGlobalParameters {
    private final double meanBias;
    private final double biasVariance;
    private final double outlierProbability;

    AlleleFractionGlobalParameters(final double meanBias,
                                   final double biasVariance,
                                   final double outlierProbability) {
        this.meanBias = meanBias;
        this.biasVariance = biasVariance;
        this.outlierProbability = outlierProbability;
    }

    double getMeanBias() {
        return meanBias;
    }

    double getBiasVariance() {
        return biasVariance;
    }

    double getOutlierProbability() {
        return outlierProbability;
    }

    //get the gamma distribution alpha parameter
    double getAlpha() {
        return meanBias * meanBias / biasVariance;
    }

    //get the gamma distribution beta parameter
    double getBeta() {
        return meanBias / biasVariance;
    }

    AlleleFractionGlobalParameters copyWithNewMeanBias(final double newMeanBias) {
        return new AlleleFractionGlobalParameters(newMeanBias, biasVariance, outlierProbability);
    }

    AlleleFractionGlobalParameters copyWithNewBiasVariance(final double newBiasVariance) {
        return new AlleleFractionGlobalParameters(meanBias, newBiasVariance, outlierProbability);
    }

    AlleleFractionGlobalParameters copyWithNewOutlierProbability(final double newOutlierProbability) {
        return new AlleleFractionGlobalParameters(meanBias, biasVariance, newOutlierProbability);
    }

    @Override
    public String toString() {
        return "AlleleFractionGlobalParameters{" +
                "meanBias=" + CopyNumberFormatsUtils.formatDouble(meanBias) +
                ", biasVariance=" + CopyNumberFormatsUtils.formatDouble(biasVariance) +
                ", outlierProbability=" + CopyNumberFormatsUtils.formatDouble(outlierProbability) +
                '}';
    }
}
