package org.broadinstitute.hellbender.engine;

import org.broadinstitute.hellbender.engine.MultiIntervalShard;
import org.broadinstitute.hellbender.engine.Shard;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.Serializable;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

/**
 * adapts a normal Shard into a MultiIntervalShard that contains only the single wrapped shard
 *
 * this is a temporary shim until we can fully adopt MultiIntervalShard in HaplotypeCallerSpark
 */
public class ShardToMultiIntervalShardAdapter<T> implements MultiIntervalShard<T>, Shard<T>, Serializable {
    private static final long serialVersionUID = 0L;
    private final Shard<T> shard;

    public ShardToMultiIntervalShardAdapter(final Shard<T> shard) {
        this.shard = Utils.nonNull(shard);
    }

    @Override
    public List<SimpleInterval> getIntervals() {
        return Collections.singletonList(shard.getInterval());
    }

    @Override
    public List<SimpleInterval> getPaddedIntervals() {
        return Collections.singletonList(shard.getPaddedInterval());
    }

    @Override
    public Iterator<T> iterator() {
        return shard.iterator();
    }

    @Override
    public SimpleInterval getInterval() {
        return shard.getInterval();
    }

    @Override
    public SimpleInterval getPaddedInterval() {
        return shard.getPaddedInterval();
    }
}
