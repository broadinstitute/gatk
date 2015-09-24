package org.broadinstitute.hellbender.engine.dataflow.datasources;

import org.broadinstitute.hellbender.engine.dataflow.datasources.ReadContextData;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.collections.IntervalsSkipListOneContig;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.Variant;

import java.io.Serializable;
import java.util.ArrayList;

/**
 *
 * Immutable holder class, just for reads.
 * (we have a separate class so we can assign it our shiny coder)
 * This class enforces no consistency guarantees. The intent
 * is for the reads to start in the interval, but that's up the callers.
 *
 * Use ReadsShardCoder instead of Java serialization for this one.
 */
public class ReadsShard implements Serializable {
    private static final long serialVersionUID = 1L;

    // the interval covered by this shard
    public final SimpleInterval interval;
    // reads that start in the shard
    public final ArrayList<GATKRead> reads;

    public ReadsShard(SimpleInterval interval) {
        this.interval = interval;
        this.reads = new ArrayList<>();
    }

    /**
     * Careful: this ctor takes ownership of the passed reads.
     * Do not modify them after this call (ideally don't even keep a reference to them).
     */
    public ReadsShard(SimpleInterval interval,final ArrayList<GATKRead> reads) {
        this.interval = interval;
        this.reads = reads;
    }

}
