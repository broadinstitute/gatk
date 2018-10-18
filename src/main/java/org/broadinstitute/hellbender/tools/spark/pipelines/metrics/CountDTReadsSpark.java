package org.broadinstitute.hellbender.tools.spark.pipelines.metrics;

import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CoverageAnalysisProgramGroup;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.filters.WellformedReadFilter;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.io.PrintStream;
import java.util.Arrays;
import java.util.List;

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
public final class CountDTReadsSpark extends GATKSparkTool {

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
    public List<ReadFilter> getDefaultReadFilters() {
        return Arrays.asList(new ReadFilterLibrary.AllowAllReadsReadFilter());
    }

    @Override
    protected void runTool(final JavaSparkContext ctx) {
        final JavaRDD<GATKRead> reads = getReads();

        final long count = reads.count();
        System.out.println("Total reads: "+count);

        final long dtReadscount = reads.filter(r -> r.hasAttribute("DT")).count();
        System.out.println("Reads marked as OpticalDuplicates: "+dtReadscount);

        final long mappeddtReadscount = reads.filter(r -> r.hasAttribute("DT")).filter(r -> !ReadUtils.isNonPrimary(r)).count();
        System.out.println("Mapped Reads marked as OpticalDuplicates: "+mappeddtReadscount);

        final long mappedReadsFirstInpair = reads.filter(r -> r.hasAttribute("DT")).filter(r -> !ReadUtils.isNonPrimary(r)).filter(r -> r.isFirstOfPair()).count();
        System.out.println("Mapped Reads first in pair marked as OpticalDuplicates: "+mappedReadsFirstInpair);

        final long readnamesRepresented = reads.filter(r -> r.hasAttribute("DT")).map(r -> r.getName()).distinct().count();
        System.out.println("Unique readnames represented: "+readnamesRepresented);

        if(out != null) {
            writeReads(ctx, out, reads.filter(r -> r.hasAttribute("DT")), getHeaderForReads());
        }
    }
}