package org.broadinstitute.hellbender.tools.copynumber.gcnv;

import org.apache.commons.math3.util.FastMath;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Map;

/**
 * A record containing a copy number posterior distribution for a single interval
 */
public class CopyNumberPosteriorRecord {

    /**
     * Posterior distribution stores in log scale
     */
    private final Map<IntegerCopyNumberState, Double> copyNumberPosteriorDistribution;

    /**
     * Tolerated error for checking equality of doubles
     */
    private final static double TOLERATED_ERROR = 1e-4;

    /**
     * @param copyNumberPosteriorDistribution copy number posterior distribution in log scale
     */
    public CopyNumberPosteriorRecord(final Map<IntegerCopyNumberState, Double> copyNumberPosteriorDistribution) {
        this.copyNumberPosteriorDistribution = Utils.nonNull(copyNumberPosteriorDistribution);
        if (MathUtils.compareDoubles(this.copyNumberPosteriorDistribution.values().stream().map(p -> FastMath.exp(p))
                .reduce((p1, p2) -> p1 + p2).get(), 1.0, TOLERATED_ERROR) != 0) {
            throw new IllegalArgumentException("Posterior probabilities for at at least one posterior record do not sum up to one");
        }
    }

    /**
     * Get the probability for a given copy number state
     */
    public double getCopyNumberPosterior(final IntegerCopyNumberState integerCopyNumberState) {
        return copyNumberPosteriorDistribution.get(integerCopyNumberState);
    }
}
