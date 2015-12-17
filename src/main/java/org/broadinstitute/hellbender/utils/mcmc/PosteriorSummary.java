package org.broadinstitute.hellbender.utils.mcmc;

/**
 * Represents the mean and standard deviation of the posterior of a univariate model parameter.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class PosteriorSummary {
    private final double mean;
    private final double standardDeviation;

    public PosteriorSummary(final double mean, final double standardDeviation) {
        this.mean = mean;
        this.standardDeviation = standardDeviation;
    }

    public PosteriorSummary(final PosteriorSummary summary) {
        this(summary.mean(), summary.standardDeviation());
    }

    public double mean() {
        return mean;
    }

    public double standardDeviation() {
        return standardDeviation;
    }
}
