package org.broadinstitute.hellbender.tools.sv;

import htsjdk.tribble.Feature;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Objects;

public final class BafEvidence implements Feature {
    final String sample;
    final String contig;
    final int position;
    final double value;

    public final static String BCI_VERSION = "1.0";

    public BafEvidence(final String sample, final String contig, final int position, final double value) {
        Utils.nonNull(sample);
        Utils.nonNull(contig);
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
        return position;
    }

    public double getValue() {
        return value;
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