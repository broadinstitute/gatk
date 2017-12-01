package org.broadinstitute.hellbender.tools.exome.segmentation;

import org.broadinstitute.hellbender.utils.param.ParamUtils;

/**
 * Created by davidben on 9/6/16.
 */
public final class AFCRHiddenState {
    private final double minorAlleleFraction;
    private final double log2CopyRatio;

    public AFCRHiddenState(double minorAlleleFraction, double log2CopyRatio) {
        this.minorAlleleFraction = ParamUtils.inRange(minorAlleleFraction, 0, 0.5, "Minor fraction must be in [0, 1/2].");
        this.log2CopyRatio = log2CopyRatio;
    }

    public double getMinorAlleleFraction() {
        return minorAlleleFraction;
    }

    public double getLog2CopyRatio() {
        return log2CopyRatio;
    }

    @Override
    public String toString() {
        return String.format("(maf=%f, log2cr=%f)", minorAlleleFraction, log2CopyRatio);
    }
}
