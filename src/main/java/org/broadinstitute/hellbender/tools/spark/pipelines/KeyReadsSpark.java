package org.broadinstitute.hellbender.tools.spark.pipelines;

import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.SparkProgramGroup;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import scala.Tuple2;

import java.io.Serializable;
import java.util.Random;

@CommandLineProgramProperties(summary = "Sorts the input BAM", oneLineSummary = "Sorts a BAM file", programGroup = SparkProgramGroup.class)

public final class KeyReadsSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;

    @Override
    public boolean requiresReads() { return true; }

    @Override
    protected void runTool(final JavaSparkContext ctx) {
        JavaRDD<GATKRead> reads = getReads();
        JavaPairRDD<String, Integer> keyed = reads.mapToPair(read -> new Tuple2<>(Utils.key(read), read.getStart()));
        JavaPairRDD<String, Iterable<Integer>> grouped = keyed.groupByKey();
        System.out.println(grouped.count());
    }
}
class Utils implements Serializable {
    private static final long serialVersionUID = 1L;

    static String key(GATKRead r) {
        try {
            if (r.isUnmapped()) {
                return r.getFragmentLength() + "," + r.getName().hashCode();
            }
            return r.getContig() + "," + ReadUtils.getStrandedUnclippedStart(r) + "," + ReadUtils.readHasMappedMate(r);
        } catch (java.lang.IllegalArgumentException e) {
            throw new RuntimeException("found the bad record: " + r.getContig() + "," + r.getStart());
        }
    }
}
