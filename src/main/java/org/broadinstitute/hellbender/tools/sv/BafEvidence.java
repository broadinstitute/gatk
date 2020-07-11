package org.broadinstitute.hellbender.tools.sv;

import htsjdk.tribble.Feature;

public final class BafEvidence implements Feature {

    final String sample;
    final String contig;
    final int position;
    final double value;

    public BafEvidence(final String sample, final String contig, final int position, final double value) {
        this.sample = sample;
        this.contig = contig;
        this.position = position;
        this.value = value;
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

    public double getValue() {
        return value;
    }
}