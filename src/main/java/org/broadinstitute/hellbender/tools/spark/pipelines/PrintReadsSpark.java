package org.broadinstitute.hellbender.tools.spark.pipelines;

import htsjdk.samtools.SAMFileHeader;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.SparkProgramGroup;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.read.GATKRead;

@CommandLineProgramProperties(summary = "Print reads from the input BAM", oneLineSummary = "PrintReads on Spark", programGroup = SparkProgramGroup.class)
public final class PrintReadsSpark extends GATKSparkTool {

    private static final long serialVersionUID = 1L;

    @Override
    public boolean requiresReads() { return true; }

    @Argument(doc = "uri for the output file: a local file path",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            optional = false)
    public String output;

    @Override
    protected void runTool(final JavaSparkContext ctx) {
        if (getHeaderForReads().getSortOrder() != SAMFileHeader.SortOrder.coordinate){
            //https://github.com/broadinstitute/hellbender/issues/929
            throw new UserException("PrintReadsSpark: Only coordinate-sorted files are currently supported");
        }

        final JavaRDD<GATKRead> reads = getReads();
        writeReads(ctx, output, reads);
    }
}
