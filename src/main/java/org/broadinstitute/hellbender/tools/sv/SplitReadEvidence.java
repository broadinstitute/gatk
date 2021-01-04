package org.broadinstitute.hellbender.tools.sv;

import htsjdk.tribble.Feature;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.codecs.SplitReadEvidenceCodec;

import java.util.Arrays;
import java.util.List;
import java.util.Objects;

public final class SplitReadEvidence implements Feature {

    final String sample;
    final String contig;
    final int position;
    final int count;
    final boolean strand;

    public SplitReadEvidence(final String sample, final String contig, final int position, final int count, final boolean strand) {
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

    public String toString() {
        final List<String> data = Arrays.asList(
                contig,
                Integer.toString(position - 1),
                strand ? SplitReadEvidenceCodec.DIRECTION_RIGHT : SplitReadEvidenceCodec.DIRECTION_LEFT,
                Integer.toString(count),
                sample
        );
        return String.join(SplitReadEvidenceCodec.COL_DELIMITER, data);
    }
}