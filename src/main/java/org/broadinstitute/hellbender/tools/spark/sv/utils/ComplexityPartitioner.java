package org.broadinstitute.hellbender.tools.spark.sv.utils;

import org.apache.spark.Partitioner;

import java.util.Arrays;

/** A Spark Partitioner that puts tasks with greater complexities into earlier partitions. */
public final class ComplexityPartitioner extends Partitioner {
    private static final long serialVersionUID = 1L;
    private final int[] partitions;

    public ComplexityPartitioner( final int[] complexities ) {
        final Integer[] tags = new Integer[complexities.length];
        for ( int idx = 0; idx != tags.length; ++idx ) {
            tags[idx] = idx;
        }
        Arrays.sort(tags, (a, b) -> Integer.compare(complexities[b], complexities[a]));
        partitions = new int[complexities.length];
        for ( int idx = 0; idx != tags.length; ++idx ) {
            partitions[tags[idx]] = idx;
        }
    }

    @Override public int numPartitions() { return partitions.length; }
    @Override public int getPartition( final Object key ) { return partitions[(Integer)key]; }
}
