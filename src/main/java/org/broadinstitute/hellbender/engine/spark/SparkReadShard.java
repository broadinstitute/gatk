package org.broadinstitute.hellbender.engine.spark;


import org.broadinstitute.hellbender.engine.Shard;
import org.broadinstitute.hellbender.engine.ShardBoundary;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.Serializable;
import java.util.Iterator;

/**
 * A simple shard implementation intended to be used for splitting reads by partition in Spark tools
 */
public final class SparkReadShard implements Shard<GATKRead>, Serializable {
    private static final long serialVersionUID = 1L;

    private final ShardBoundary boundaries;
    private final Iterable<GATKRead> reads;

    public SparkReadShard(final ShardBoundary boundaries, final Iterable<GATKRead> reads){
        this.boundaries = Utils.nonNull(boundaries);
        this.reads = Utils.nonNull(reads);
    }

    @Override
    public SimpleInterval getInterval() {
        return boundaries.getInterval();
    }

    @Override
    public SimpleInterval getPaddedInterval() {
        return boundaries.getPaddedInterval();
    }

    @Override
    public Iterator<GATKRead> iterator() {
        return reads.iterator();
    }
}
