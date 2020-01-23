package org.broadinstitute.hellbender.tools.walkers.mutect.clustering;

import org.apache.commons.math3.special.Gamma;
import org.broadinstitute.hellbender.tools.walkers.mutect.SomaticLikelihoodsEngine;
import org.broadinstitute.hellbender.tools.walkers.readorientation.BetaDistributionShape;
import org.broadinstitute.hellbender.tools.walkers.validation.basicshortmutpileup.BetaBinomialDistribution;

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
    public double correctedLogLikelihood(final Datum datum) {
        return correctedLogLikelihood(datum, betaDistributionShape);
    }

    @Override
    public double logLikelihood(final int totalCount, final int altCount) {
        return new BetaBinomialDistribution(null, betaDistributionShape.getAlpha(), betaDistributionShape.getBeta(), totalCount).logProbability(altCount);
    }

    public static double correctedLogLikelihood(final Datum datum, final BetaDistributionShape betaDistributionShape) {
        final int altCount = datum.getAltCount();
        final int refCount = datum.getTotalCount() - altCount;
        return datum.getTumorLogOdds() + logOddsCorrection(BetaDistributionShape.FLAT_BETA, betaDistributionShape, altCount, refCount);
    }

    @Override
    public void learn(final List<Datum> data, final double[] responsibilities) {
        double alpha = betaDistributionShape.getAlpha();
        double beta = betaDistributionShape.getBeta();

        for (int epoch = 0; epoch < NUM_EPOCHS; epoch++) {
            for (int n = 0; n < data.size(); n++) {
                final Datum datum = data.get(n);
                final int alt = datum.getAltCount();
                final int ref = datum.getTotalCount() - alt;

                final double digammaOfTotalPlusAlphaPlusBeta = Gamma.digamma(datum.getTotalCount() + alpha + beta);
                final double digammaOfAlphaPlusBeta = Gamma.digamma(alpha + beta);
                final double alphaGradient = Gamma.digamma(alpha + alt) - digammaOfTotalPlusAlphaPlusBeta - Gamma.digamma(alpha) + digammaOfAlphaPlusBeta;
                final double betaGradient = Gamma.digamma(beta + ref) - digammaOfTotalPlusAlphaPlusBeta - Gamma.digamma(beta) + digammaOfAlphaPlusBeta;

                alpha = Math.max(alpha + RATE * alphaGradient * responsibilities[n], 1.0);
                beta = Math.max(beta + RATE * betaGradient * responsibilities[n], 0.5);
            }
        }

        betaDistributionShape = new BetaDistributionShape(alpha, beta);

    }

    private static double logOddsCorrection(final BetaDistributionShape originalBeta, final BetaDistributionShape newBeta, final int altCount, final int refCount) {
        return g(newBeta.getAlpha(), newBeta.getBeta()) - g(newBeta.getAlpha() + altCount, newBeta.getBeta() + refCount)
                - g(originalBeta.getAlpha(), originalBeta.getBeta()) + g(originalBeta.getAlpha() + altCount, originalBeta.getBeta() + refCount);
    }

    private static double g(final double... omega) {
        return SomaticLikelihoodsEngine.logDirichletNormalization(omega);
    }

    @Override
    public String toString() {
        return String.format("alpha = %.2f, beta = %.2f", betaDistributionShape.getAlpha(), betaDistributionShape.getBeta());
    }
}
