package org.broadinstitute.hellbender.engine.spark;

import org.apache.spark.RangePartitioner;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.engine.spark.datasources.ReadsSparkSink;
import org.broadinstitute.hellbender.utils.read.GATKRead;

public class TestRangePartitionAndWrite extends GATKSparkTool {

    @Argument
    boolean repartition = false;
    @Override
    protected void runTool(JavaSparkContext ctx) {
        final JavaRDD<GATKRead> reads = getReads();
        JavaRDD<GATKRead> output;
        if(repartition) {
            final JavaPairRDD<String, GATKRead> keyedReads = reads.keyBy(GATKRead::getName);
            final RangePartitioner<Object, Object> partitioner = new RangePartitioner<>();
            final JavaPairRDD<String, GATKRead> partitioned = keyedReads.partitionBy(partitioner);
            output = partitioned.map(kv -> kv._2);
        } else {
            output = reads;
        }
        writeReads(ctx, "out.bam", output );
    }
}
