package org.broadinstitute.hellbender.utils.variant;

import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.io.Serializable;

/**
 * SkeletonVariant is a minimal implementation of the Variant interface.
 */
public class SkeletonVariant implements Variant, Serializable {
    private static final long serialVersionUID = 1L;

    final private SimpleInterval interval;
    final private boolean snp;
    final private boolean indel;

    public SkeletonVariant(SimpleInterval interval, boolean isSNP, boolean isIndel) {
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
        return String.format("SkeletonVariant -- interval(%s:%d-%d), snp(%b), indel(%b)",
                interval.getContig(), interval.getStart(), interval.getEnd(), snp, indel);
    }
}
