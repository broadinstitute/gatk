package org.broadinstitute.hellbender.tools.sv;

import org.broadinstitute.hellbender.utils.Utils;

import java.util.Objects;
import java.util.Set;

/** Biallelic-frequency of a sample at some locus. */
public final class BafEvidence implements SVFeature {
    private final String sample;
    private final String contig;
    private final int position;
    private final double value;

    public final static String BCI_VERSION = "1.0";

    public BafEvidence( final String sample, final String contig,
                        final int position, final double value ) {
        Utils.nonNull(sample);
        Utils.nonNull(contig);
        this.sample = sample;
        this.contig = contig;
        this.position = position;
        this.value = value;
    }

    // value-altering constructor
    public BafEvidence( final BafEvidence that, final double value ) {
        this.sample = that.sample;
        this.contig = that.contig;
        this.position = that.position;
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
    public BafEvidence extractSamples( final Set<String> sampleList, final Object header ) {
        return sampleList.contains(sample) ? this : null;
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

    @Override public String toString() {
        return contig + "\t" + position + "\t" + sample + "\t" + value;
    }
}
