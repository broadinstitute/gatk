package org.broadinstitute.hellbender.tools.spark.sv;

import htsjdk.samtools.SAMSequenceDictionary;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.SparkProgramGroup;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import scala.Tuple2;

import java.util.ArrayList;
import java.util.List;

@CommandLineProgramProperties(summary="Computes binned coverage and average mapq",
        oneLineSummary="Computes binned coverage and average mapq",
        programGroup = SparkProgramGroup.class)
public class ComputeBinnedCoverageAndMapqSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;

    @Argument(doc = "URL of the output path", shortName = "outputPath",
            fullName = "outputPath", optional = false)
    private String outputPath;

    @Override
    public boolean requiresReads() {
        return true;
    }

    @Override
    public boolean requiresReference() {
        return true;
    }

    @Override
    protected void runTool(final JavaSparkContext ctx) {
        final SAMSequenceDictionary referenceSequenceDictionary = getReferenceSequenceDictionary();
        final JavaPairRDD<Tuple2<Integer, Integer>, Tuple2<Integer, Integer>> readCountAndTotalQualByBin = getReads().
                filter(r -> !r.isSecondaryAlignment() && !r.isSupplementaryAlignment() && !r.failsVendorQualityCheck() && !r.isUnmapped()).
                flatMapToPair(gatkRead -> {
                    Integer contig = referenceSequenceDictionary.getSequenceIndex(gatkRead.getContig());
                    List<Tuple2<Tuple2<Integer, Integer>, Tuple2<Integer, Integer>>> results = new ArrayList<>();
                    int binSize = 100;
                    for (int i = gatkRead.getStart(); i < gatkRead.getStart() + gatkRead.getLength(); i = i + binSize) {
                        results.add(new Tuple2<>(new Tuple2<>(contig, i - i % binSize), new Tuple2<>(1, gatkRead.getMappingQuality())));
                    }
                    return results;
                }).reduceByKey((p1, p2) -> new Tuple2<>(p1._1() + p2._1(), p1._2() + p2._2()));

        readCountAndTotalQualByBin.saveAsTextFile(outputPath);

    }
}
