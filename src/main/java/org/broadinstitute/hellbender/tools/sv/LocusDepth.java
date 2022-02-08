package org.broadinstitute.hellbender.tools.sv;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.util.Locatable;

import java.util.*;

/** The read depth of each base call for a sample at some locus. */
@VisibleForTesting
public final class LocusDepth implements SVFeature {
    private final String contig;
    private final int position;
    private final String sample;
    private final byte refCall; // index into nucleotideValues
    private final int[] depths;
    public final static String BCI_VERSION = "1.0";

    public LocusDepth( final Locatable loc, final String sample, final byte refCall ) {
        this.contig = loc.getContig();
        this.position = loc.getStart();
        this.sample = sample;
        this.refCall = refCall;
        this.depths = new int[4];
    }

    public LocusDepth( final String contig, final int position,
                       final String sample, final byte refCall,
                       final int aDepth, final int cDepth, final int gDepth, final int tDepth ) {
        this.contig = contig;
        this.position = position;
        this.sample = sample;
        this.refCall = refCall;
        this.depths = new int[4];
        depths[0] = aDepth;
        depths[1] = cDepth;
        depths[2] = gDepth;
        depths[3] = tDepth;
    }

    public void observe( final int idx ) {
        depths[idx] += 1;
    }

    @Override
    public String getContig() {
        return contig;
    }

    @Override
    public int getEnd() {
        return position;
    }

    @Override
    public int getStart() {
        return position;
    }

    public String getSample() { return sample; }
    public char getRefCall() {
        return (char)refCall;
    }
    public int getADepth() { return depths[0]; }
    public int getCDepth() { return depths[1]; }
    public int getGDepth() { return depths[2]; }
    public int getTDepth() { return depths[3]; }

    @Override
    public LocusDepth extractSamples( final Set<String> sampleNames, final Object header ) {
        return sampleNames.contains(sample) ? this : null;
    }

    @Override
    public boolean equals( final Object obj ) {
        if ( !(obj instanceof LocusDepth) ) return false;
        return equals((LocusDepth)obj);
    }

    public boolean equals( final LocusDepth that ) {
        return this.position == that.position && this.refCall == that.refCall &&
                this.depths[0] == that.depths[0] && this.depths[1] == that.depths[1] &&
                this.depths[2] == that.depths[2] && this.depths[3] == that.depths[3] &&
                Objects.equals(this.contig, that.contig) &&
                Objects.equals(this.sample, that.sample);
    }

    @Override
    public int hashCode() {
        return Objects.hash(contig, position, sample, refCall, Arrays.hashCode(depths));
    }

    @Override
    public String toString() {
        return contig + "\t" + position + "\t" + sample + "\t" + (char)refCall + "\t" +
                depths[0] + "\t" + depths[1] + "\t" + depths[2] + "\t" + depths[3];
    }
}
