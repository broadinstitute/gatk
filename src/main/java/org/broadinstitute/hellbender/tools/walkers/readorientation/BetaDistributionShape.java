package org.broadinstitute.hellbender.tools.walkers.readorientation;

import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

public class BetaDistributionShape {
    public static final BetaDistributionShape FLAT_BETA = new BetaDistributionShape(1,1);

    private double alpha;
    private double beta;

    public BetaDistributionShape(final double alpha, final double beta){
        ParamUtils.isPositive(alpha, "alpha must be greater than 0 but got " + alpha);
        Utils.validateArg(beta > 0, "beta must be greater than 0 but got " + beta);

        this.alpha = alpha;
        this.beta = beta;
    }

    public double getAlpha() {
        return alpha;
    }

    public double getBeta() {
        return beta;
    }
}
