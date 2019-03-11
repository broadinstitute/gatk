package org.broadinstitute.hellbender.tools.walkers.mutect.clustering;

import org.apache.commons.math3.special.Gamma;
import org.broadinstitute.hellbender.tools.walkers.mutect.SomaticLikelihoodsEngine;
import org.broadinstitute.hellbender.tools.walkers.readorientation.BetaDistributionShape;
import org.broadinstitute.hellbender.tools.walkers.validation.basicshortmutpileup.BetaBinomialDistribution;
import org.broadinstitute.hellbender.utils.MathUtils;

import java.util.List;

public class BetaBinomialCluster implements AlleleFractionCluster {

    private static final double RATE = 0.01;
    private static final double MAX_RATE = 0.1;
    private static final int NUM_EPOCHS = 10;

    BetaDistributionShape betaDistributionShape;

    public BetaBinomialCluster(final BetaDistributionShape betaDistributionShape) {
        this.betaDistributionShape = betaDistributionShape;
    }

    @Override
    public double log10Likelihood(final Datum datum) {
        return log10Likelihood(datum, betaDistributionShape);
    }

    @Override
    public double log10Likelihood(final int totalCount, final int altCount) {
        return MathUtils.LOG10_OF_E * new BetaBinomialDistribution(null, betaDistributionShape.getAlpha(), betaDistributionShape.getBeta(), totalCount).logProbability(altCount);
    }

    public static double log10Likelihood(final Datum datum, final BetaDistributionShape betaDistributionShape) {
        final int altCount = datum.getAltCount();
        final int refCount = datum.getTotalCount() - altCount;
        return datum.getTumorLog10Odds() + log10OddsCorrection(BetaDistributionShape.FLAT_BETA, betaDistributionShape, altCount, refCount);
    }

    @Override
    public void learn(final List<Datum> data) {
        double alpha = betaDistributionShape.getAlpha();
        double beta = betaDistributionShape.getBeta();

        for (int epoch = 0; epoch < NUM_EPOCHS; epoch++) {
            for (final Datum datum : data) {
                final int alt = datum.getAltCount();
                final int ref = datum.getTotalCount() - alt;

                final double digammaOfTotalPlusAlphaPlusBeta = Gamma.digamma(datum.getTotalCount() + alpha + beta);
                final double digammaOfAlphaPlusBeta = Gamma.digamma(alpha + beta);
                final double alphaGradient = Gamma.digamma(alpha + alt) - digammaOfTotalPlusAlphaPlusBeta - Gamma.digamma(alpha) + digammaOfAlphaPlusBeta;
                final double betaGradient = Gamma.digamma(beta + ref) - digammaOfTotalPlusAlphaPlusBeta - Gamma.digamma(beta) + digammaOfAlphaPlusBeta;

                alpha = Math.max(alpha + RATE * alphaGradient, 0.5);
                beta = Math.max(beta + RATE * betaGradient, 0.5);
            }
        }

        betaDistributionShape = new BetaDistributionShape(alpha, beta);

    }

    private static double log10OddsCorrection(final BetaDistributionShape originalBeta, final BetaDistributionShape newBeta, final int altCount, final int refCount) {
        return g(newBeta.getAlpha(), newBeta.getBeta()) - g(newBeta.getAlpha() + altCount, newBeta.getBeta() + refCount)
                - g(originalBeta.getAlpha(), originalBeta.getBeta()) + g(originalBeta.getAlpha() + altCount, originalBeta.getBeta() + refCount);
    }

    private static double g(final double... omega) {
        return SomaticLikelihoodsEngine.log10DirichletNormalization(omega);
    }

    @Override
    public String toString() {
        return String.format("alpha = %.2f, beta = %.2f", betaDistributionShape.getAlpha(), betaDistributionShape.getBeta());
    }
}
