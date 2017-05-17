package org.broadinstitute.hellbender.utils.locusiterator;

import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.downsampling.DownsampleType;
import org.broadinstitute.hellbender.utils.downsampling.DownsamplingMethod;

import java.io.Serializable;

/**
 * Simple wrapper about the information LIBS needs about downsampling
 */
public final class LIBSDownsamplingInfo implements Serializable {
    private static final long serialVersionUID = 1L;

    private final boolean performDownsampling;
    private final int toCoverage;

    /**
     * @param performDownsampling whether to downsample
     * @param toCoverage what coverage to downsample to (or -1 if no downsampling)
     */
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

    public static LIBSDownsamplingInfo toDownsamplingInfo(final DownsamplingMethod downsamplingMethod) {
        final boolean performDownsampling = downsamplingMethod != null &&
                downsamplingMethod.type == DownsampleType.BY_SAMPLE &&
                downsamplingMethod.toCoverage != null;

        final int toCoverage = performDownsampling ? downsamplingMethod.toCoverage : 0;

        return new LIBSDownsamplingInfo(performDownsampling, toCoverage);
    }
}
