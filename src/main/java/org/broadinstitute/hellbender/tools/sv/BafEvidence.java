package org.broadinstitute.hellbender.tools.sv;

import htsjdk.tribble.Feature;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.codecs.BafEvidenceCodec;

import java.util.Arrays;
import java.util.List;
import java.util.Objects;

/**
 * B-Allele frequency (BAF) evidence container
 */
public final class BafEvidence implements Feature {

    final String sample;
    final String contig;
    final int position;
    final double value;

    public BafEvidence(final String sample, final String contig, final int position, final double value) {
        this.sample = Utils.nonNull(sample);
        this.contig = Utils.nonNull(contig);
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

    @Override
    public String toString() {
        final List<String> data = Arrays.asList(
                contig,
                Integer.toString(position - 1),
                Double.toString(value),
                sample
        );
        return String.join(BafEvidenceCodec.COL_DELIMITER, data);
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof BafEvidence)) return false;
        BafEvidence that = (BafEvidence) o;
        return position == that.position &&
                Double.compare(that.value, value) == 0 &&
                sample.equals(that.sample) &&
                contig.equals(that.contig);
    }

    @Override
    public int hashCode() {
        return Objects.hash(sample, contig, position, value);
    }
}