package org.broadinstitute.hellbender.engine;

import java.io.Serializable;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Iterator;

/**
 * A {@link Shard} backed by a {@link ShardBoundary} and a collection of records.
 */
public final class ShardBoundaryShard<T> implements Shard<T>, Serializable {
    private static final long serialVersionUID = 1L;

    private final ShardBoundary shardBoundary;
    private final Iterable<T> locatables;

    /**
     * Create a new {@link ShardBoundaryShard} from the given {@link ShardBoundary} and records.
     * @param shardBoundary the boundary defining the shard
     * @param locatables the records overlapping  the shard
     */
    public ShardBoundaryShard(final ShardBoundary shardBoundary, final Iterable<T> locatables) {
        Utils.nonNull(this.shardBoundary = shardBoundary);
        Utils.nonNull(this.locatables = locatables);
    }

    @Override
    public SimpleInterval getInterval() {
        return shardBoundary.getInterval();
    }

    @Override
    public SimpleInterval getPaddedInterval() {
        return shardBoundary.getPaddedInterval();
    }

    @Override
    public ShardBoundary getBoundary() {
        return shardBoundary;
    }

    @Override
    public Iterator<T> iterator() {
        return locatables.iterator();
    }
}
