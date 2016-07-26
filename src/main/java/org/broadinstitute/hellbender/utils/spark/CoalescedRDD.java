package org.broadinstitute.hellbender.utils.spark;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Iterators;
import com.google.common.primitives.Ints;

import java.util.List;

import org.apache.spark.Dependency;
import org.apache.spark.NarrowDependency;
import org.apache.spark.Partition;
import org.apache.spark.TaskContext;
import org.apache.spark.rdd.PartitionGroup;
import org.apache.spark.rdd.RDD;
import scala.collection.Iterator;
import scala.collection.JavaConversions;
import scala.collection.Seq;
import scala.reflect.ClassTag;
import scala.reflect.ClassTag$;

class CoalescedRDD<T> extends RDD<T> {

    private static final long serialVersionUID = 1L;

    private transient RDD<T> prev;
    private final int maxPartitions;
    private final PartitionCoalescer partitionCoalescer;
    private Class<T> cls;

    public CoalescedRDD(RDD<T> prev, int maxPartitions, PartitionCoalescer
            partitionCoalescer, Class<T> cls) {
        super(prev.context(), null, ClassTag$.MODULE$.apply(cls));

        this.prev = prev;
        this.maxPartitions = maxPartitions;
        this.partitionCoalescer = partitionCoalescer;
        this.cls = cls;
    }

    @Override
    public Partition[] getPartitions() {
        PartitionGroup[] partitionGroups = partitionCoalescer.coalesce(maxPartitions, prev);
        Partition[] partitions = new Partition[partitionGroups.length];
        for (int i = 0; i < partitions.length; i++) {
            PartitionGroup pg = partitionGroups[i];
            List<Partition> partitionList = JavaConversions.asJavaList(partitionGroups[i].arr());
            int[] ids = partitionList.stream().mapToInt(Partition::index).toArray();
            partitions[i] = new CoalescedRDDPartition(i, prev, ids, pg.prefLoc());
        }
        return partitions;
    }

    @Override
    public Iterator<T> compute(Partition split, TaskContext context) {
        java.util.Iterator<java.util.Iterator<T>> iterators =
                ((CoalescedRDDPartition) split).getParents().stream().map(p -> {
                            ClassTag<T> tag = ClassTag$.MODULE$.apply(cls);
                            RDD<T> objectRDD = firstParent(tag);
                            return JavaConversions.asJavaIterator(objectRDD.iterator(p, context));
                        }
                ).iterator();
        return JavaConversions.asScalaIterator(Iterators.concat(iterators));
    }

    @Override
    @SuppressWarnings({"unchecked", "rawtypes"})
    public Seq<Dependency<?>> getDependencies() {
        return JavaConversions.asScalaBuffer(ImmutableList.of(
                new NarrowDependency(prev) {
                    private static final long serialVersionUID = 1L;
                    @Override
                    public Seq<Integer> getParents(int partitionId) {
                        List<Integer> i = Ints.asList(
                                ((CoalescedRDDPartition) partitions()[partitionId]).getParentsIndices());
                        return JavaConversions.asScalaBuffer(i);
                    }
                }
        ));
    }

    @Override
    public void clearDependencies() {
        super.clearDependencies();
        prev = null;
    }

    @Override
    public Seq<String> getPreferredLocations(Partition partition) {
        return ((CoalescedRDDPartition) partition).getPreferredLocation().toList();
    }
}
