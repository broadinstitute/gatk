package org.broadinstitute.hellbender.utils.variant;

import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.io.Serializable;

/**
 * MinimalVariant is a minimal implementation of the {@link GATKVariant} interface.
 */
public class MinimalVariant implements GATKVariant, Serializable {
    private static final long serialVersionUID = 1L;

    private final SimpleInterval interval;
    private final boolean snp;
    private final boolean indel;

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        MinimalVariant that = (MinimalVariant) o;

        if (isSnp() != that.isSnp()) return false;
        if (isIndel() != that.isIndel()) return false;
        return !(interval != null ? !interval.equals(that.interval) : that.interval != null);

    }

    @Override
    public int hashCode() {
        int result = interval != null ? interval.hashCode() : 0;
        result = 31 * result + (isSnp() ? 1 : 0);
        result = 31 * result + (isIndel() ? 1 : 0);
        return result;
    }

    public MinimalVariant(SimpleInterval interval, boolean isSNP, boolean isIndel) {
        this.interval = interval;
        this.snp = isSNP;
        this.indel = isIndel;
    }

    @Override
    public String getContig() { return interval.getContig(); }
    @Override
    public int getStart() { return interval.getStart(); }
    @Override
    public int getEnd() { return interval.getEnd(); }
    @Override
    public boolean isSnp() { return snp; }
    @Override
    public boolean isIndel() { return indel; }

    @Override
    public String toString() {
        return String.format("MinimalVariant -- interval(%s:%d-%d), snp(%b), indel(%b)",
                interval.getContig(), interval.getStart(), interval.getEnd(), snp, indel);
    }
}
