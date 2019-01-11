package org.broadinstitute.hellbender.engine.spark;

import org.broadinstitute.barclay.argparser.Argument;

import java.io.Serializable;

public class AssemblyRegionReadShardArgumentCollection implements Serializable {
    private static final long serialVersionUID = 1L;

    public static final int DEFAULT_READSHARD_SIZE = 5000;
    public static final int DEFAULT_READSHARD_PADDING_SIZE = 100;

    @Argument(fullName="readShardSize", shortName="readShardSize", doc = "Maximum size of each read shard, in bases. For good performance, this should be much larger than the maximum assembly region size.", optional = true)
    public int readShardSize = DEFAULT_READSHARD_SIZE;

    @Argument(fullName="readShardPadding", shortName="readShardPadding", doc = "Each read shard has this many bases of extra context on each side. Read shards must have as much or more padding than assembly regions.", optional = true)
    public int readShardPadding = DEFAULT_READSHARD_PADDING_SIZE;
}
