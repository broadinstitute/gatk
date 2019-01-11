package org.broadinstitute.hellbender.tools.spark.pipelines;

import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import picard.cmdline.programgroups.DiagnosticsAndQCProgramGroup;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.tools.FlagStat.FlagStatus;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.PrintStream;

/**
 * Spark tool to accumulate flag statistics given a BAM file, e.g. total number of reads with QC failure flag set, number of
 * duplicates, percentage mapped etc.
 *
 * <p>
 * This tool can be run without explicitly specifying Spark options. That is to say,
 * the given example command without Spark options will run locally.
 * See <a href ="https://software.broadinstitute.org/gatk/documentation/article?id=10060">Tutorial#10060</a>
 * for an example of how to set up and run a Spark tool on a cloud Spark cluster.
 * </p>
 *
 * <h3>Input</h3>
 * <ul>
 *     <li>A BAM file containing aligned read data</li>
 * </ul>
 *
 * <h3>Output</h3>
 * <ul>
 *     <li>Accumulated flag statistics</li>
 * </ul>
 *
 * <h3>Example Usage</h3>
 * <pre>
 *   gatk FlagStatSpark \
 *     -I input.bam \
 *     -O statistics.txt
 * </pre>
 */
@BetaFeature
@DocumentedFeature
@CommandLineProgramProperties(
        summary = "Spark tool to accumulate flag statistics given a BAM file, e.g. total number of reads with QC failure flag set," +
                "number of duplicates, percentage mapped etc.",
        oneLineSummary = "Spark tool to accumulate flag statistics",
        programGroup = DiagnosticsAndQCProgramGroup.class
)
public final class FlagStatSpark extends GATKSparkTool {

    private static final long serialVersionUID = 1L;

    @Override
    public boolean requiresReads() { return true; }

    @Argument(
            doc = "uri for the output file: a local file path",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            optional = true
    )
    public String out;

    @Override
    protected void runTool(final JavaSparkContext ctx) {
        final JavaRDD<GATKRead> reads = getReads();

        final FlagStatus result = reads.aggregate(new FlagStatus(), FlagStatus::add, FlagStatus::merge);
        System.out.println(result);

        if(out != null ) {
            try ( final PrintStream ps = new PrintStream(BucketUtils.createFile(out)) ) {
                ps.print(result);
            }
        }
    }
}
