package org.broadinstitute.hellbender.tools.sv;

import org.broadinstitute.hellbender.utils.Utils;

import java.util.Objects;
import java.util.Set;

/**
 * Documents evidence of reads (of some sample at some locus) that align well to reference for
 * some portion of the read, and fails to align for another portion of the read.
 * Strand actually refers to whether the fails-to-align bit is upstream (at the beginning of the read)
 * or downstream (at the end of the read) relative to the part that aligns.  I think that
 * strand==true, encoded as "right", means that the non-aligned part is at the end of the read.
 * Unless maybe it means the opposite.
 */
public final class SplitReadEvidence implements SVFeature {

    private final String sample;
    private final String contig;
    private final int position;
    private final int count;
    private final boolean strand;

    public final static String BCI_VERSION = "1.0";

    public SplitReadEvidence( final String sample, final String contig, final int position,
                              final int count, final boolean strand ) {
        Utils.nonNull(sample);
        Utils.nonNull(contig);
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
        return position;
    }

    public boolean getStrand() {
        return strand;
    }

    public int getCount() {
        return count;
    }

    @Override
    public SplitReadEvidence extractSamples( final Set<String> sampleNames, final Object header ) {
        return sampleNames.contains(sample) ? this : null;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof SplitReadEvidence)) return false;
        SplitReadEvidence that = (SplitReadEvidence) o;
        return position == that.position &&
                count == that.count &&
                strand == that.strand &&
                sample.equals(that.sample) &&
                contig.equals(that.contig);
    }

    @Override
    public int hashCode() {
        return Objects.hash(sample, contig, position, count, strand);
    }

    @Override public String toString() {
        return contig + "\t" + position + "\t" + sample + "\t" + count + "\t" + strand;
    }
}
