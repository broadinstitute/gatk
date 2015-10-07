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
     * Use a shuffle join strategy, where both sides of join are shuffled across the workers.
     */
    SHUFFLE
}
