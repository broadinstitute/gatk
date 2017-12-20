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
 * Calculate the overall number of bases SAM/BAM/CRAM file
 *
 * <h3>Input</h3>
 * <ul>
 *     <li>A single BAM file</li>
 * </ul>
 *
 * <h3>Output</h3>
 * <ul>
 *     <li>A text file containing number of bases</li>
 * </ul>
 *
 * <h3>Example</h3>
 *
 * <h4>Output base count to file</h4>
 * <pre>
 *   gatk CountBasesSpark \
 *     -I input_reads.bam \
 *     -O base_count.txt
 * </pre>
 *
 * <h4>Print base count</h4>
 * <pre>
 *     gatk CountBasesSpark \
 *       -I input_reads.bam
 * </pre>
 */
@CommandLineProgramProperties(
        summary = "Counts bases in the input SAM/BAM",
        oneLineSummary = "Counts bases in the input SAM/BAM",
        programGroup = CoverageAnalysisProgramGroup.class
)
@DocumentedFeature
@BetaFeature
public final class CountBasesSpark extends GATKSparkTool {

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

        final long count = reads.map(r -> (long)r.getLength()).reduce(Long::sum);
        System.out.println(count);

        if( out != null) {
            try ( final PrintStream ps = new PrintStream(BucketUtils.createFile(out)) ) {
                ps.print(count);
            }
        }
    }
}
