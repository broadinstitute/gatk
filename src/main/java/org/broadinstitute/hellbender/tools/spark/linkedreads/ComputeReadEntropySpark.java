package org.broadinstitute.hellbender.tools.spark.linkedreads;

import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.SparkProgramGroup;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.tools.spark.sv.SVKmer;
import org.broadinstitute.hellbender.tools.spark.sv.SVKmerizer;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.HashMap;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;

@CommandLineProgramProperties(
        summary = "Computes read entropy",
        oneLineSummary = "ComputeReadEntropy on Spark",
        programGroup = SparkProgramGroup.class
)
public class ComputeReadEntropySpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;

    @Argument(doc = "uri for the output file",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            optional = true)
    public String out;

    @Override
    public boolean requiresReads() {
        return true;
    }

    @Override
    protected void runTool(final JavaSparkContext ctx) {

        final JavaRDD<GATKRead> reads = getReads();

        reads.mapToDouble(r -> computeEntropy(r.getBasesString())).stats();

    }

    public static double computeEntropy(String sequence) {
        final Map<SVKmer, Long> frequencies = SVKmerizer.stream(sequence, 3).collect(Collectors.groupingBy(k -> k, Collectors.counting()));
        final int numKmers = sequence.length() - 2;
        final Double entropy = -1 * frequencies.entrySet().stream().collect(Collectors.summingDouble(e -> {
            final double p = (double) e.getValue() / numKmers;
            return p * Math.log(p) / Math.log(2.0d);
        }));
        return entropy;
    }
}