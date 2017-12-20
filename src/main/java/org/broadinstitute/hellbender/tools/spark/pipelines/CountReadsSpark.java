package org.broadinstitute.hellbender.tools.spark.pipelines;

import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CoverageAnalysisProgramGroup;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.PrintStream;

/**
 * Calculate the overall number of reads in a SAM/BAM file
 *
 * <h3>Input</h3>
 * <ul>
 *     <li>A single BAM file</li>
 * </ul>
 *
 * <h3>Output</h3>
 * <ul>
 *     <li>A text file containing number of reads</li>
 * </ul>
 *
 * <h3>Example</h3>
 *
 * <h4>Output number of reads to file</h4>
 * <pre>
 *   gatk CountReadsSpark \
 *     -I input_reads.bam \
 *     -O read_count.txt
 * </pre>
 *
 * <h4>Print read count</h4>
 * <pre>
 *   gatk CountReadsSpark \
 *     -I input_reads.bam
 * </pre>
 */
@DocumentedFeature
@BetaFeature
@CommandLineProgramProperties(
        summary = "Counts reads in the input SAM/BAM",
        oneLineSummary = "Counts reads in the input SAM/BAM",
        programGroup = CoverageAnalysisProgramGroup.class
)
public final class CountReadsSpark extends GATKSparkTool {

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

        final long count = reads.count();
        System.out.println(count);

        if(out != null) {
            try (final PrintStream ps = new PrintStream(BucketUtils.createFile(out))) {
                ps.print(count);
            }
        }
    }
}
