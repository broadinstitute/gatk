package org.broadinstitute.hellbender.tools.sv;

import org.broadinstitute.hellbender.utils.codecs.BafEvidenceCodec;

import java.util.Arrays;
import java.util.List;

public final class BafEvidence extends SVEvidence {

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

    public String encode() {
        final List<String> data = Arrays.asList(
                contig,
                Integer.toString(position - 1),
                Double.toString(value),
                sample
        );
        return String.join(BafEvidenceCodec.COL_DELIMITER, data);
    }
}