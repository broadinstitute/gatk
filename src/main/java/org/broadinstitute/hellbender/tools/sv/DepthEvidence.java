package org.broadinstitute.hellbender.tools.sv;

import htsjdk.tribble.Feature;

public final class DepthEvidence implements Feature {

    final String contig;
    final int start;
    final int end;
    final int[] counts;

    public DepthEvidence(final String contig, int start, final int end, final int[] counts) {
        this.contig = contig;
        this.start = start;
        this.end = end;
        this.counts = counts;
    }

    @Override
    public String getContig() {
        return contig;
    }

    @Override
    public int getStart() {
        return start;
    }

    @Override
    public int getEnd() {
        return end;
    }

    public int[] getCounts() { return counts; }
}