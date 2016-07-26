package org.broadinstitute.hellbender.utils.spark;

import org.apache.spark.rdd.PartitionGroup;
import org.apache.spark.rdd.RDD;

interface PartitionCoalescer {
    PartitionGroup[] coalesce(int maxPartitions, RDD<?> parent);
}
