package org.broadinstitute.hellbender.tools.sv;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;

/** The read depth of each base call for a sample at some locus. */
@VisibleForTesting
public final class LocusDepth implements SVFeature {
    private final String contig;
    private final int position;
    private final String sample;
    // next three fields index "ACGT"
    private final int refIndex;
    private final int altIndex;
    private final int[] depths;

    public final static String BCI_VERSION = "1.0";

    public LocusDepth( final Locatable loc, final String sample, final int refIndex, final int altIndex ) {
        Utils.validateArg(refIndex >= 0 && refIndex <= 3, "refIndex must be between 0 and 3");
        Utils.validateArg(altIndex >= 0 && altIndex <= 3, "altIndex must be between 0 and 3");
        Utils.validateArg(refIndex != altIndex, "refIndex and altIndex must be different");
        this.contig = loc.getContig();
        this.position = loc.getStart();
        this.sample = sample;
        this.refIndex = refIndex;
        this.altIndex = altIndex;
        this.depths = new int[4];
    }

    public LocusDepth( final String contig, final int position,
                       final String sample, final int refIndex, final int altIndex,
                       final int aDepth, final int cDepth, final int gDepth, final int tDepth ) {
        Utils.validateArg(refIndex >= 0 && refIndex <= 3, "refIndex must be between 0 and 3");
        Utils.validateArg(altIndex >= 0 && altIndex <= 3, "altIndex must be between 0 and 3");
        Utils.validateArg(refIndex != altIndex, "refIndex and altIndex must be different");
        this.contig = contig;
        this.position = position;
        this.sample = sample;
        this.refIndex = refIndex;
        this.altIndex = altIndex;
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
    // index into "ACGT"
    public int getRefIndex() {
        return refIndex;
    }
    // index into "ACGT"
    public int getAltIndex() {
        return altIndex;
    }
    public int getTotalDepth() { return depths[0] + depths[1] + depths[2] + depths[3]; }
    // idx is index into "ACGT"
    public int getDepth( final int idx ) { return depths[idx]; }
    public int getRefDepth() { return depths[refIndex]; }
    public int getAltDepth() { return depths[altIndex]; }

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
        return this.position == that.position && this.refIndex == that.refIndex &&
                this.depths[0] == that.depths[0] && this.depths[1] == that.depths[1] &&
                this.depths[2] == that.depths[2] && this.depths[3] == that.depths[3] &&
                Objects.equals(this.contig, that.contig) &&
                Objects.equals(this.sample, that.sample);
    }

    @Override
    public int hashCode() {
        return Objects.hash(contig, position, sample, refIndex, Arrays.hashCode(depths));
    }

    @Override
    public String toString() {
        return contig + "\t" + position + "\t" + sample + "\t" + (char)refIndex + "\t" +
                depths[0] + "\t" + depths[1] + "\t" + depths[2] + "\t" + depths[3];
    }
}
