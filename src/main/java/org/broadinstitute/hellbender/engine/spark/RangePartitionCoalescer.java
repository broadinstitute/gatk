package org.broadinstitute.hellbender.engine.spark;

import org.apache.spark.Partition;
import org.apache.spark.rdd.PartitionCoalescer;
import org.apache.spark.rdd.PartitionGroup;
import org.apache.spark.rdd.RDD;
import scala.collection.JavaConversions;
import scala.collection.Seq;

import java.io.Serializable;
import java.util.Arrays;
import java.util.List;

/**
 * A {@link PartitionCoalescer} that allows a range of partitions to be coalesced into groups.
 */
class RangePartitionCoalescer implements PartitionCoalescer, Serializable, scala.Serializable {

    private static final long serialVersionUID = 1L;

    private List<Integer> maxEndPartitionIndexes;

    /**
     * @param maxEndPartitionIndexes the indexes of the end of each coalesced partition, so that
     *                               coalesced partition {@code i} in the coalesced partitions is made up of partitions
     *                               from index {@code i} to {@code maxEndPartitionIndexes.get(i)} (inclusive)
     */
    public RangePartitionCoalescer(List<Integer> maxEndPartitionIndexes) {
        this.maxEndPartitionIndexes = maxEndPartitionIndexes;
    }

    @Override
    public PartitionGroup[] coalesce(int maxPartitions, RDD<?> parent) {
        if (maxPartitions != parent.getNumPartitions()) {
            throw new IllegalArgumentException("Cannot use " + getClass().getSimpleName() +
                    " with a different number of partitions to the parent RDD.");
        }
        List<Partition> partitions = Arrays.asList(parent.getPartitions());
        PartitionGroup[] groups = new PartitionGroup[partitions.size()];

        for (int i = 0; i < partitions.size(); i++) {
            Seq<String> preferredLocations = parent.getPreferredLocations(partitions.get(i));
            scala.Option<String> preferredLocation = scala.Option.apply
                    (preferredLocations.isEmpty() ? null : preferredLocations.apply(0));
            PartitionGroup group = new PartitionGroup(preferredLocation);
            List<Partition> partitionsInGroup =
                    partitions.subList(i, maxEndPartitionIndexes.get(i) + 1);
            group.partitions().append(JavaConversions.asScalaBuffer(partitionsInGroup));
            groups[i] = group;
        }
        return groups;
    }
}
