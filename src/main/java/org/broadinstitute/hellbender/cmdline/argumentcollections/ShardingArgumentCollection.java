package org.broadinstitute.hellbender.cmdline.argumentcollections;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineException;

import java.io.Serializable;

/**
 * Argument collection for shard intervals.
 *
 * Note: it does not include any specific AssemblyRegionWalker parameter, just parameters to generate {@link org.broadinstitute.hellbender.engine.Shard}.
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public abstract class ShardingArgumentCollection implements Serializable {
    private static final long serialVersionUID = 1L;

    /**
     * If readShardSize is set to this value, we will not shard the user's intervals. Instead,
     * we'll create one shard per interval (or one shard per contig, if intervals are not explicitly specified)
     */
    public static final int NO_INTERVAL_SHARDING = -1;

    public abstract int getShardSize();

    public abstract int getShardPadding();

    public abstract int getShardStep();

    /**
     * Validates the arguments.
     * @throws CommandLineException.BadArgumentValue if they are not valid
     */
    public final void validateArguments() {
        if ( getShardSize() <= 0 && getShardSize() != NO_INTERVAL_SHARDING) {
            throw new CommandLineException.BadArgumentValue("read shard size must be > 0 or " + NO_INTERVAL_SHARDING + " for no sharding");
        }

        if ( getShardStep() <= 0 ) {
            throw new CommandLineException.BadArgumentValue("shard step must be > 0");
        }

        if ( getShardPadding() < 0 ) {
            throw new CommandLineException.BadArgumentValue("shard padding must be >= 0");
        }
    }

}
