package org.broadinstitute.hellbender.utils.variant.writers;

public final class ReblockingOptions {

    private boolean dropLowQuals = false;
    private boolean allowMissingHomRefData = false;
    private double rgqThreshold = 0.0;

    public ReblockingOptions() {}

    public ReblockingOptions(final boolean dropLowQuals, final boolean allowMissingHomRefData, final double rgqThreshold) {
        this.dropLowQuals = dropLowQuals;
        this.allowMissingHomRefData = allowMissingHomRefData;
        this.rgqThreshold = rgqThreshold;
    }

    public boolean getDropLowQualsOpt() {
        return dropLowQuals;
    }

    public boolean getAllowMissingHomRefDataOpt() {
        return allowMissingHomRefData;
    }

    public double getRgqThresholdOpt() {
        return rgqThreshold;
    }
}
