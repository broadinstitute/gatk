package org.broadinstitute.hellbender.engine.spark;

import org.apache.spark.rdd.PartitionGroup;
import org.apache.spark.rdd.RDD;

/**
 * A PartitionCoalescer defines how to coalesce the partitions of a given RDD.
 * A backport of Spark's PartitionCoalescer from Spark 2.0 to work with Spark 1.6.
 */
interface PartitionCoalescer {
    /**
     * Coalesce the partitions of the given RDD.
     *
     * @param maxPartitions the maximum number of partitions to have after coalescing
     * @param parent the parent RDD whose partitions to coalesce
     * @return an array of [[PartitionGroup]]s, where each element is itself an array of
     * [[Partition]]s and represents a partition after coalescing is performed.
     */
    PartitionGroup[] coalesce(int maxPartitions, RDD<?> parent);
}
