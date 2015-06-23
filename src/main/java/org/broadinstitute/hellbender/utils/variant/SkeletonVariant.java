package org.broadinstitute.hellbender.utils.variant;

import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.io.Serializable;
import java.util.UUID;

/**
 * SkeletonVariant is a minimal implementation of the Variant interface.
 */
public class SkeletonVariant implements Variant, Serializable {
    private static final long serialVersionUID = 1L;

    final private SimpleInterval interval;
    final private boolean snp;
    final private boolean indel;
    final private UUID uuid;

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        SkeletonVariant that = (SkeletonVariant) o;

        if (isSnp() != that.isSnp()) return false;
        if (isIndel() != that.isIndel()) return false;
        if (!interval.equals(that.interval)) return false;
        return uuid.equals(that.uuid);

    }

    @Override
    public int hashCode() {
        int result = interval.hashCode();
        result = 31 * result + (isSnp() ? 1 : 0);
        result = 31 * result + (isIndel() ? 1 : 0);
        result = 31 * result + uuid.hashCode();
        return result;
    }

    public SkeletonVariant(SimpleInterval interval, boolean isSNP, boolean isIndel, UUID uuid) {
        this.interval = interval;
        this.snp = isSNP;
        this.indel = isIndel;
        this.uuid = uuid;
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
    public UUID getUUID() {
        return uuid;
    }

    @Override
    public String toString() {
        return String.format("SkeletonVariant -- interval(%s:%d-%d), snp(%b), indel(%b)",
                interval.getContig(), interval.getStart(), interval.getEnd(), snp, indel);
    }
}
