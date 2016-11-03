package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.Serializable;
import java.util.Objects;

/**
 * Holds the bounds of a {@link Shard}, both with and without padding
 */
public class ShardBoundary implements Locatable, Serializable {
    private static final long serialVersionUID = 1L;

    private final SimpleInterval interval;
    private final SimpleInterval paddedInterval;

    /**
     * Create a new ShardBoundary from the given intervals
     *
     * @param interval the interval covered by the shard
     * @param paddedInterval the interval covered by the shard's padding, must contain the shard interval
     */
    public ShardBoundary(final SimpleInterval interval, final SimpleInterval paddedInterval) {
        Utils.nonNull(interval);
        Utils.nonNull(paddedInterval);
        Utils.validateArg(paddedInterval.contains(interval), "interval must be contained within paddedInterval");
        this.interval = interval;
        this.paddedInterval = paddedInterval;
    }


    @Override
    public String getContig() {
        return interval.getContig();
    }

    /**
     * @return start of the shard boundary without padding
     */
    @Override
    public int getStart() {
        return interval.getStart();
    }

    /**
     * @return end of the shard boundary without padding
     */
    @Override
    public int getEnd() {
        return interval.getEnd();
    }

    /**
     * @return the interval this boundary covers without padding
     */
    public SimpleInterval getInterval() {
        return interval;
    }

    /**
     * @return the interval this boundary covers including padding
     */
    public SimpleInterval getPaddedInterval() {
        return paddedInterval;
    }

    @Override
    public boolean equals(Object o) {
        if ( this == o ) {
            return true;
        }
        if ( o == null || getClass() != o.getClass() ) {
            return false;
        }

        final ShardBoundary key = (ShardBoundary) o;

        if ( !Objects.equals(interval, key.interval) ) {
            return false;
        }
        return Objects.equals(paddedInterval, key.paddedInterval);

    }

    @Override
    public int hashCode() {
        return Objects.hash(interval, paddedInterval);
    }
}
