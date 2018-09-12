package org.broadinstitute.hellbender.tools.funcotator;

/**
 * Simple struct container class for the 5'/3' flank settings. Helps prevent the two values from
 * getting swapped accidentally when passing them around.
 */
public final class FlankSettings {

    public final int fivePrimeFlankSize;
    public final int threePrimeFlankSize;

    public FlankSettings(final int fivePrimeFlankSize, final int threePrimeFlankSize) {
        this.fivePrimeFlankSize = fivePrimeFlankSize;
        this.threePrimeFlankSize = threePrimeFlankSize;
    }
}
