package org.broadinstitute.hellbender.tools.spark.sv.discovery.readdepth;

import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;

public final class SVCopyRatio {
    private final int contig;
    private final int start;
    private final int end;
    private final float log2CopyRatio;

    public SVCopyRatio(final int contig, final int start, final int end, final float log2CopyRatio) {
        this.contig = contig;
        this.start = start;
        this.end = end;
        this.log2CopyRatio = log2CopyRatio;
    }

    public SVCopyRatio(final SVInterval interval, final float log2CopyRatio) {
        this.contig = interval.getContig();
        this.start = interval.getStart();
        this.end = interval.getEnd();
        this.log2CopyRatio = log2CopyRatio;
    }

    public int getContig() {
        return contig;
    }

    public int getStart() {
        return start;
    }

    public int getEnd() {
        return end;
    }

    public float getLog2CopyRatio() {
        return log2CopyRatio;
    }

    public SVInterval getInterval() {
        return new SVInterval(contig, start, end);
    }
}
