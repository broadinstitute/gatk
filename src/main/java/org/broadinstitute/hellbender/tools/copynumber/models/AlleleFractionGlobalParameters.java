package org.broadinstitute.hellbender.tools.copynumber.models;

/**
 * Encapsulates the global parameters of the allele fraction model: the mean and variance of the common prior on
 * allelic biases and the outlier probability.
 *
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 */
final class AlleleFractionGlobalParameters {
    static final String DOUBLE_FORMAT = MultidimensionalModeller.DOUBLE_FORMAT;

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
        return String.format("AlleleFractionGlobalParameters{" +
                "meanBias=" + DOUBLE_FORMAT +
                ", biasVariance=" + DOUBLE_FORMAT +
                ", outlierProbability=" + DOUBLE_FORMAT +
                '}', meanBias, biasVariance, outlierProbability);
    }
}
