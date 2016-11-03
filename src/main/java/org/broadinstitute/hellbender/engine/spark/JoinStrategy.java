package org.broadinstitute.hellbender.engine.spark;

/**
 * Possible join strategies when using Spark
 */
public enum JoinStrategy {
    /**
     * Use a broadcast join strategy, where one side of the join is collected into memory and broadcast to all workers.
     */
    BROADCAST,

    /**
     * Use an overlaps partitioner strategy, where one side of the join is sharded in partitions and the other side is broadcast.
     */
    OVERLAPS_PARTITIONER,

    /**
     * Use a shuffle join strategy, where both sides of join are shuffled across the workers.
     */
    SHUFFLE
}
