package org.broadinstitute.hellbender.engine;

import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.Objects;

/**
 * A cut-down version of {@link AssemblyRegion} that doesn't store reads, used in the strict implementation of
 * {@link org.broadinstitute.hellbender.engine.spark.FindAssemblyRegionsSpark}.
 */
public class ReadlessAssemblyRegion extends ShardBoundary {
    private static final long serialVersionUID = 1L;

    private final boolean isActive;

    public ReadlessAssemblyRegion(final AssemblyRegion assemblyRegion) {
        super(assemblyRegion.getSpan(), assemblyRegion.getPaddedSpan());
        this.isActive = assemblyRegion.isActive();
    }

    private ReadlessAssemblyRegion(final SimpleInterval activeRegionLoc, SimpleInterval paddedSpan, final boolean isActive, final boolean padded) {
        super(activeRegionLoc, paddedSpan, padded);
        this.isActive = isActive;
    }

    public boolean isActive() {
        return isActive;
    }

    @Override
    public ShardBoundary paddedShardBoundary() {
        return padded ? this : new ReadlessAssemblyRegion(interval, paddedSpan, isActive, true);
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        if (!super.equals(o)) return false;
        ReadlessAssemblyRegion that = (ReadlessAssemblyRegion) o;
        return isActive == that.isActive;
    }

    @Override
    public int hashCode() {
        return Objects.hash(super.hashCode(), isActive);
    }

    @Override
    public String toString() {
        return "ReadlessAssemblyRegion{" +
                "isActive=" + isActive +
                ", interval=" + interval +
                ", paddedSpan=" + paddedSpan +
                ", padded=" + padded +
                '}';
    }
}
