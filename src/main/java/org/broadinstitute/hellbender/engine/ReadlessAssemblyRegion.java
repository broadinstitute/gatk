package org.broadinstitute.hellbender.engine;

import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.Objects;

// TODO: see if it's possible to make the serialized representation more compact for the Spark shuffle

/**
 * A cut-down version of {@link AssemblyRegion} that doesn't store reads, used in
 * {@link org.broadinstitute.hellbender.engine.spark.NewAssemblyRegionWalkerSpark}.
 */
public class ReadlessAssemblyRegion extends ShardBoundary {
    private static final long serialVersionUID = 1L;

    private final boolean isActive;

    public ReadlessAssemblyRegion(final AssemblyRegion assemblyRegion) {
        super(assemblyRegion.getSpan(), assemblyRegion.getExtendedSpan());
        this.isActive = assemblyRegion.isActive();
    }

    private ReadlessAssemblyRegion(final SimpleInterval activeRegionLoc, SimpleInterval extendedLoc, final boolean isActive, final boolean padded) {
        super(activeRegionLoc, extendedLoc, padded);
        this.isActive = isActive;
    }

    public boolean isActive() {
        return isActive;
    }

    @Override
    public ShardBoundary paddedShardBoundary() {
        return padded ? this : new ReadlessAssemblyRegion(interval, paddedInterval, isActive, true);
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
                ", paddedInterval=" + paddedInterval +
                ", padded=" + padded +
                '}';
    }
}
