package org.broadinstitute.hellbender.utils.spark;

import org.apache.spark.Partition;
import org.apache.spark.rdd.RDD;
import scala.Option;

import java.io.IOException;
import java.io.ObjectOutputStream;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

class CoalescedRDDPartition implements Partition {

    private static final long serialVersionUID = 1L;

    private int index;
    private transient RDD<?> rdd;
    private int[] parentsIndices;
    private transient Option<String> preferredLocation;
    private List<Partition> parents;

    public CoalescedRDDPartition(int index, RDD<?> rdd, int[] parentsIndices,
                                 Option<String> preferredLocation) {
        this.index = index;
        this.rdd = rdd;
        this.parentsIndices = parentsIndices;
        this.preferredLocation = preferredLocation;
        parents = Arrays.stream(parentsIndices)
                .mapToObj(i -> rdd.partitions()[i])
                .collect(Collectors.toList());
    }

    @Override
    public int index() {
        return index;
    }

    public int[] getParentsIndices() {
        return parentsIndices;
    }

    public Option<String> getPreferredLocation() {
        return preferredLocation;
    }

    public List<Partition> getParents() {
        return parents;
    }

    private void writeObject(ObjectOutputStream oos) throws IOException {
        // Update the reference to parent partition at the time of task serialization
        parents = Arrays.stream(parentsIndices)
                .mapToObj(i -> rdd.partitions()[i])
                .collect(Collectors.toList());
        oos.defaultWriteObject();
    }
}
