package org.broadinstitute.hellbender.tools.spark.pipelines;

import htsjdk.samtools.SAMFileHeader;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.SparkProgramGroup;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.engine.spark.datasources.ReadsSparkSink;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadsWriteFormat;

import java.io.IOException;

@CommandLineProgramProperties(summary = "Print reads from the input BAM", oneLineSummary = "Print reads from the input BAM", programGroup = SparkProgramGroup.class)
public final class PrintReadsSpark extends GATKSparkTool {

    public static final String SHARDED_OUTPUT_LONG_ARG = "shardedOutput";
    public static final String SHARDED_OUTPUT_SHORT_ARG = "shardedOutput";

    private static final long serialVersionUID = 1L;

    @Override
    public boolean requiresReads() { return true; }

    @Argument(doc = "uri for the output file: a local file path",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            optional = false)
    public String output;

    @Argument(doc = "If specified, shard the output bam", shortName = SHARDED_OUTPUT_SHORT_ARG, fullName = SHARDED_OUTPUT_LONG_ARG, optional = true)
    public boolean shardedOutput = false;

    @Override
    protected void runTool(final JavaSparkContext ctx) {
        if (getHeaderForReads().getSortOrder() != SAMFileHeader.SortOrder.coordinate){
            //https://github.com/broadinstitute/hellbender/issues/929
            throw new UserException("PrintReadsSpark: Only coordinate-sorted files are currently supported");
        }

        final JavaRDD<GATKRead> reads = getReads();

        try {
            ReadsSparkSink.writeReads(ctx, output, reads, getHeaderForReads(), shardedOutput ? ReadsWriteFormat.SHARDED : ReadsWriteFormat.SINGLE);
        } catch (final IOException e) {
            throw new GATKException("unable to write bam: " + e);
        }
    }
}