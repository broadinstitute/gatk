package org.broadinstitute.hellbender.tools.walkers.mutect.clustering;

import org.broadinstitute.hellbender.tools.walkers.readorientation.BetaDistributionShape;
import org.broadinstitute.hellbender.tools.walkers.validation.basicshortmutpileup.BetaBinomialDistribution;
import org.broadinstitute.hellbender.utils.MathUtils;

import java.util.List;
import java.util.stream.IntStream;

public class BinomialCluster implements AlleleFractionCluster {

    private static final double STD_DEV_OVER_MEAN = 0.01;

    private BetaDistributionShape betaDistributionShape;

    public BinomialCluster(final double mean) {
        betaDistributionShape = getFuzzyBinomial(mean);
    }

    @Override
    public double correctedLogLikelihood(final Datum datum) {
        return BetaBinomialCluster.correctedLogLikelihood(datum, betaDistributionShape);
    }

    @Override
    public double logLikelihood(final int totalCount, final int altCount) {
        return new BetaBinomialDistribution(null, betaDistributionShape.getAlpha(), betaDistributionShape.getBeta(), totalCount).logProbability(altCount);
    }

    @Override
    public void learn(final List<Datum> data, final double[] responsibilities) {
        final double altCount = IntStream.range(0, data.size()).mapToDouble(n -> data.get(n).getAltCount()*responsibilities[n]).sum() + 0.0001;
        final double totalCount = IntStream.range(0, data.size()).mapToDouble(n -> data.get(n).getTotalCount()*responsibilities[n]).sum() + 0.0001;
        betaDistributionShape = getFuzzyBinomial( altCount / totalCount);

    }

    private static BetaDistributionShape getFuzzyBinomial(final double unboundedMean) {
        final double mean = Math.min(unboundedMean, 1 - STD_DEV_OVER_MEAN);
        final double alphaPlusBeta = ((1 - mean) / (mean * MathUtils.square(STD_DEV_OVER_MEAN))) - 1;
        final double alpha = mean * alphaPlusBeta;
        final double beta = alphaPlusBeta - alpha;
        return new BetaDistributionShape(alpha, beta);
    }

    @Override
    public String toString() {
        return String.format("mean = %.3f", betaDistributionShape.getAlpha() / (betaDistributionShape.getAlpha() + betaDistributionShape.getBeta()));
    }
}
