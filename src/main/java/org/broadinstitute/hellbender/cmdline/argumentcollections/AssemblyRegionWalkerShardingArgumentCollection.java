package org.broadinstitute.hellbender.cmdline.argumentcollections;

import org.broadinstitute.barclay.argparser.Argument;

/**
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public class AssemblyRegionWalkerShardingArgumentCollection extends ShardingArgumentCollection {
    private static final long serialVersionUID = 1L;


    @Argument(fullName="readShardSize", shortName="readShardSize", doc = "Maximum size of each read shard, in bases. For good performance, this should be much larger than the maximum assembly region size.", optional = true)
    protected int readShardSize;

    @Argument(fullName="readShardPadding", shortName="readShardPadding", doc = "Each read shard has this many bases of extra context on each side. Read shards must have as much or more padding than assembly regions.", optional = true)
    protected int readShardPadding;

    public AssemblyRegionWalkerShardingArgumentCollection(final int defaultReadShardSize, final int defaultReadShardPadding) {
        this.readShardSize = defaultReadShardSize;
        this.readShardPadding = defaultReadShardPadding;
    }

    @Override
    public int getShardSize() {
        return readShardSize;
    }

    @Override
    public int getShardPadding() {
        return readShardPadding;
    }

    @Override
    public int getShardStep() {
        return readShardSize;
    }
}
