package org.broadinstitute.hellbender.engine.dataflow.datasources;

import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.collections.IntervalsSkipListOneContig;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.Variant;

import java.io.Serializable;
import java.util.ArrayList;

/**
 * Immutable holder class.
 * This class enforces no consistency guarantees. The intent of course
 * is for the reads and variants to be in the interval, but that's up the callers.
 */
public class ContextShard implements Serializable {
    private static final long serialVersionUID = 1L;

    // the interval covered by this shard
    public final SimpleInterval interval;
    // variants that overlap with the shard.
    public final IntervalsSkipListOneContig<Variant> variants;
    // reads that start in the shard
    public final ArrayList<GATKRead> reads;
    // variants and reference for the particular read at the same index as this element.
    public final ArrayList<ReadContextData> readContext;

    public ContextShard(SimpleInterval interval) {
        this.interval = interval;
        this.variants = null;
        this.reads = new ArrayList<>();
        this.readContext = new ArrayList<>();
    }

    /**
     * Careful: this ctor takes ownership of the passed reads and ReadContextData array.
     * Do not modify them after this call (ideally don't even keep a reference to them).
     */
    private ContextShard(SimpleInterval interval, IntervalsSkipListOneContig<Variant> variants, final ArrayList<GATKRead> reads, final ArrayList<ReadContextData> readContext) {
        this.interval = interval;
        this.variants = variants;
        this.reads = reads;
        this.readContext = readContext;
    }

    /**
     * create a new shard, keeping only the variants that overlap
     * with the new interval. Reads etc are unchanged.
     */
    public ContextShard split(SimpleInterval newInterval) {
        final IntervalsSkipListOneContig<Variant> newVariants;
        if (null==variants) {
            newVariants = null;
        } else {
            newVariants = new IntervalsSkipListOneContig<>( variants.getOverlapping(newInterval) );
        }
        return new ContextShard(newInterval, newVariants, reads, readContext);
    }

    /**
     * creates a new shard, adding the specified variants.
     * Note that readContext is unchanged (including the variants it may refer to).
     */
    public ContextShard withVariants(ArrayList<Variant> newVariants) {
        return new ContextShard(this.interval, new IntervalsSkipListOneContig<>(newVariants), reads, readContext);
    }

    /**
     * creates a new shard, adding the specified reads.
     * Careful: this call takes ownership of the passed read array. So you can't modify it after this call.
     */
    public ContextShard withReads(ArrayList<GATKRead> newReads) {
        return new ContextShard(this.interval, this.variants, newReads, readContext);
    }

    /**
     * creates a new shard, adding the specified read context and *removing the variants*.
     * Careful: this call takes ownership of the passed ReadContextData array. So you can't modify it after this call.
     */
    public ContextShard withReadContext(ArrayList<ReadContextData> newReadContext) {
        return new ContextShard(interval, null, reads, newReadContext);
    }

    /**
     * Returns the variants that overlap the query interval, in start-position order.
     */
    public ArrayList<Variant> variantsOverlapping(SimpleInterval interval) {
        return variants.getOverlapping(interval);
    }

}
