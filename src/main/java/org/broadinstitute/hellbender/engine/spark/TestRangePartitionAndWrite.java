package org.broadinstitute.hellbender.engine.spark;

import org.apache.spark.RangePartitioner;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.engine.spark.datasources.ReadsSparkSink;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import scala.math.Ordering;
import scala.math.Ordering$;
import scala.reflect.ClassTag;
import scala.reflect.ClassTag$;
import scala.reflect.ClassTag$class;

import java.util.Comparator;

public class TestRangePartitionAndWrite extends GATKSparkTool {

    private static final long serialVersionUID = 1L;
    @Argument
    boolean repartition = false;
    @Override
    protected void runTool(JavaSparkContext ctx) {
        final JavaRDD<GATKRead> reads = getReads();
        JavaRDD<GATKRead> output;
        if(repartition) {
            final JavaPairRDD<String, GATKRead> keyedReads = reads.keyBy(GATKRead::getName);

            final Ordering<String> stringOrdering = Ordering$.MODULE$.comparatorToOrdering(
                    Comparator.<String>naturalOrder());
            final ClassTag<String> classtag = ClassTag$.MODULE$.apply(Stringg.class);
            final RangePartitioner<String, GATKRead> partitioner = new RangePartitioner<String, GATKRead>(
                    reads.getNumPartitions(), keyedReads.rdd(), true,
                    stringOrdering,
                    classtag);
            final JavaPairRDD<String, GATKRead> partitioned = keyedReads.partitionBy(partitioner);
            output = partitioned.map(kv -> kv._2);
        } else {
            output = reads;
        }
        writeReads(ctx, "out.bam", output );
    }
}
