package org.broadinstitute.hellbender.tools.spark.pipelines;

import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.SparkProgramGroup;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.FlagStat.FlagStatus;
import org.broadinstitute.hellbender.utils.dataflow.BucketUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;

@CommandLineProgramProperties(summary ="runs FlagStat on Spark", oneLineSummary = "FlagStat", programGroup = SparkProgramGroup.class)
public final class FlagStatSpark extends GATKSparkTool {

    private static final long serialVersionUID = 1L;

    @Override
    public boolean requiresReads() { return true; }

    @Argument(doc = "uri for the output file: a local file path",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            optional = true)
    public String out;

    @Override
    protected void runTool(final JavaSparkContext ctx) {
        final JavaRDD<GATKRead> reads = getReads();

        final FlagStatus result = reads.aggregate(new FlagStatus(), FlagStatus::add, FlagStatus::merge);
        System.out.println(result);

        final File file = new File(out);
        try(final OutputStream outputStream = BucketUtils.createFile(file.getPath(), null);
            final PrintStream ps = new PrintStream(outputStream)) {
            ps.print(result);
        } catch(final IOException e){
            throw new UserException.CouldNotCreateOutputFile(file, e);
        }
    }
}
