package org.broadinstitute.hellbender.tools.sv;


import htsjdk.tribble.Feature;

public final class DiscordantPairEvidence implements Feature {

    final String sample;
    final String startContig;
    final String endContig;
    final int start;
    final int end;
    final boolean startStrand;
    final boolean endStrand;

    public DiscordantPairEvidence(final String sample, final String startContig, final int start, final boolean startStrand,
                                  final String endContig, final int end, final boolean endStrand) {
        this.sample = sample;
        this.startContig = startContig;
        this.start = start;
        this.startStrand = startStrand;
        this.endContig = endContig;
        this.end = end;
        this.endStrand = endStrand;
    }

    public String getSample() {
        return sample;
    }

    // For purposes of indexing, we will return the start position
    @Override
    public String getContig() {
        return startContig;
    }

    @Override
    public int getStart() {
        return start;
    }

    public boolean getStartStrand() {
        return startStrand;
    }

    public String getEndContig() {
        return endContig;
    }

    // Note this returns the end position on the start chromosome
    @Override
    public int getEnd() {
        return start + 1;
    }

    public int getEndPosition() {
        return end;
    }

    public boolean getEndStrand() {
        return endStrand;
    }
}