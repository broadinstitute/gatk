package org.broadinstitute.hellbender.utils.locusiterator;

import org.broadinstitute.hellbender.utils.Utils;

/**
 * Simple wrapper about the information LIBS needs about downsampling
 */
public final class LIBSDownsamplingInfo {
    private final boolean performDownsampling;
    private final int toCoverage;

    public LIBSDownsamplingInfo(final boolean performDownsampling, final int toCoverage) {
        Utils.validateArg(toCoverage >= -1, "toCoverage must be at least -1 (special value) but was " + toCoverage);
        this.performDownsampling = performDownsampling;
        this.toCoverage = toCoverage;
    }

    public boolean isPerformDownsampling() {
        return performDownsampling;
    }

    public int getToCoverage() {
        return toCoverage;
    }
}
