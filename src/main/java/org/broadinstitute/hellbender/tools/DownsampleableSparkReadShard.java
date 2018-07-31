package org.broadinstitute.hellbender.tools;

import java.io.Serializable;
import java.util.Iterator;
import org.broadinstitute.hellbender.engine.Shard;
import org.broadinstitute.hellbender.engine.ShardBoundary;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.downsampling.ReadsDownsampler;
import org.broadinstitute.hellbender.utils.downsampling.ReadsDownsamplingIterator;
import org.broadinstitute.hellbender.utils.read.GATKRead;

/**
 * A simple shard implementation intended to be used for splitting reads by partition in Spark tools
 */
public final class DownsampleableSparkReadShard implements Shard<GATKRead>, Serializable {
    private static final long serialVersionUID = 1L;

    private final ShardBoundary boundaries;
    private final Iterable<GATKRead> reads;
    private final ReadsDownsampler downsampler;

    /**
     * @param boundaries the boundary defining the shard
     * @param reads the records overlapping  the shard
     * @param downsampler the downsampler to use (may be null, which signifies that no downsampling is to be performed)
     */
    public DownsampleableSparkReadShard(final ShardBoundary boundaries, final Iterable<GATKRead> reads, final ReadsDownsampler downsampler) {
        this.boundaries = Utils.nonNull(boundaries);
        this.reads = Utils.nonNull(reads);
        this.downsampler = downsampler;
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
        Iterator<GATKRead> readsIterator = reads.iterator();

        if ( downsampler != null ) {
            readsIterator = new ReadsDownsamplingIterator(readsIterator, downsampler);
        }

        return readsIterator;
    }
}
