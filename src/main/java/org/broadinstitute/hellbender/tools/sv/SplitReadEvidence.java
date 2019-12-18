package org.broadinstitute.hellbender.tools.sv;

import htsjdk.tribble.Feature;
import org.broadinstitute.hellbender.utils.codecs.SplitReadEvidenceCodec;

import java.util.Arrays;
import java.util.List;

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

    @Override
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