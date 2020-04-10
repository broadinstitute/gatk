package org.broadinstitute.hellbender.tools.sv;

import htsjdk.tribble.Feature;

public final class SplitReadEvidence implements Feature {

    final String sample;
    final String contig;
    final int position;
    final int count;
    final boolean strand;

    public SplitReadEvidence(final String sample, final String contig, final int position, final int count, final boolean strand) {
        this.sample = sample;
        this.contig = contig;
        this.position = position;
        this.count = count;
        this.strand = strand;
    }

    public String getSample() {
        return sample;
    }

    @Override
    public String getContig() {
        return contig;
    }

    @Override
    public int getStart() {
        return position;
    }

    @Override
    public int getEnd() {
        return position + 1;
    }

    public boolean getStrand() {
        return strand;
    }

    public int getCount() {
        return count;
    }
}