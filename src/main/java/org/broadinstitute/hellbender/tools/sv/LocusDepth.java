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
    private final int[] depths; // four slots for each call, ACGT

    public final static String BCI_VERSION = "1.0";

    public LocusDepth( final Locatable loc, final String sample ) {
        this(loc.getContig(), loc.getStart(), sample, 0, 0, 0, 0);
    }

    public LocusDepth( final String contig, final int position, final String sample,
                       final int aDepth, final int cDepth, final int gDepth, final int tDepth ) {
        Utils.nonNull(contig, "contig may not be null");
        Utils.validateArg(position > 0, "starting position must be positive");
        Utils.nonNull(sample, "sample must not be null");
        Utils.validateArg(aDepth >= 0, "depth must be non-negative");
        Utils.validateArg(cDepth >= 0, "depth must be non-negative");
        Utils.validateArg(gDepth >= 0, "depth must be non-negative");
        Utils.validateArg(tDepth >= 0, "depth must be non-negative");
        this.contig = contig;
        this.position = position;
        this.sample = sample;
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

    public int getTotalDepth() { return depths[0] + depths[1] + depths[2] + depths[3]; }

    // idx is index into "ACGT"
    public int getDepth( final int idx ) { return depths[idx]; }

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
        return this.position == that.position &&
                this.depths[0] == that.depths[0] && this.depths[1] == that.depths[1] &&
                this.depths[2] == that.depths[2] && this.depths[3] == that.depths[3] &&
                Objects.equals(this.contig, that.contig) &&
                Objects.equals(this.sample, that.sample);
    }

    @Override
    public int hashCode() {
        return Objects.hash(contig, position, sample, Arrays.hashCode(depths));
    }

    @Override
    public String toString() {
        return sample + "@" + contig + ":" + position +
                "[" + depths[0] + "," + depths[1] + "," + depths[2] + "," + depths[3] + "]";
    }
}
