package org.broadinstitute.hellbender.engine;

import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.collections.IntervalsSkipListOneContig;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVariant;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

/**
 * Immutable storage class.
 * This class enforces no consistency guarantees. The intent
 * is for the reads to start in the interval and for ReadContext to have
 * all the variants and reference bases that overlap those reads, but that's up the callers.
 */
public class ContextShard implements Serializable {
    private static final long serialVersionUID = 1L;

    // the interval covered by this shard
    public final SimpleInterval interval;
    // variants that overlap with the shard.
    public final IntervalsSkipListOneContig<GATKVariant> variants;
    // reads that start in the shard
    public final List<GATKRead> reads;
    // variants and reference for the particular read at the same index as this element.
    public final List<ReadContextData> readContext;

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
    private ContextShard(SimpleInterval interval, IntervalsSkipListOneContig<GATKVariant> variants, final List<GATKRead> reads, final List<ReadContextData> readContext) {
        this.interval = interval;
        this.variants = variants;
        this.reads = reads;
        this.readContext = readContext;
    }

    /**
     * create a new shard, keeping only the variants that overlap
     * with the new interval. Reads, readContext, and the variants in readContext are unchanged.
     */
    public ContextShard split(SimpleInterval newInterval) {
        final IntervalsSkipListOneContig<GATKVariant> newVariants;
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
    public ContextShard withVariants(List<GATKVariant> newVariants) {
        return new ContextShard(this.interval, new IntervalsSkipListOneContig<>(newVariants), reads, readContext);
    }

    /**
     * creates a new shard, adding the specified reads.
     * Careful: this call takes ownership of the passed read array. So you can't modify it after this call.
     */
    public ContextShard withReads(List<GATKRead> newReads) {
        return new ContextShard(this.interval, this.variants, newReads, readContext);
    }

    /**
     * creates a new shard, adding the specified read context and *removing the variants*.
     * Careful: this call takes ownership of the passed ReadContextData array. So you can't modify it after this call.
     */
    public ContextShard withReadContext(List<ReadContextData> newReadContext) {
        return new ContextShard(interval, null, reads, newReadContext);
    }

    /**
     * Returns the variants that overlap the query interval, in start-position order.
     */
    public List<GATKVariant> variantsOverlapping(SimpleInterval interval) {
        return variants.getOverlapping(interval);
    }

}
