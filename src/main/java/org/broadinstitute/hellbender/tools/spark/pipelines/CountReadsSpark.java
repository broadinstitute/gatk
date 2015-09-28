package org.broadinstitute.hellbender.tools.spark.pipelines;

import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.ArgumentCollection;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.ReadInputArgumentCollection;
import org.broadinstitute.hellbender.cmdline.argumentcollections.RequiredReadInputArgumentCollection;
import org.broadinstitute.hellbender.cmdline.programgroups.SparkProgramGroup;
import org.broadinstitute.hellbender.engine.spark.SparkCommandLineProgram;
import org.broadinstitute.hellbender.engine.spark.datasources.ReadsSparkSource;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.dataflow.BucketUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;

@CommandLineProgramProperties(summary = "Counts reads in the input BAM", oneLineSummary = "Counts reads in a BAM file", programGroup = SparkProgramGroup.class)
public final class CountReadsSpark extends SparkCommandLineProgram {

    private static final long serialVersionUID = 1L;

    @ArgumentCollection
    public ReadInputArgumentCollection readArguments= new RequiredReadInputArgumentCollection();;

    @Argument(doc = "uri for the output file: a local file path",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            optional = true)
    public String out;

    @Override
    protected void runPipeline(final JavaSparkContext ctx) {
        if ( readArguments.getReadFilesNames().size() != 1 ) {
            throw new UserException("This tool only accepts a single bam/sam/cram as input");
        }
        final String bam = readArguments.getReadFilesNames().get(0);
        final ReadsSparkSource readSource = new ReadsSparkSource(ctx);
        final JavaRDD<GATKRead> reads = readSource.getParallelReads(bam, null);
        final long count = reads.count();
        System.out.println(count);
        if (out != null){
            final File file = new File(out);
            try(final OutputStream outputStream = BucketUtils.createFile(file.getPath(), null);
                final PrintStream ps = new PrintStream(outputStream)) {
                ps.print(count);
            } catch(final IOException e){
                throw new UserException.CouldNotCreateOutputFile(file, e);
            }
        }
    }

    @Override
    protected String getProgramName() {
        return getClass().getSimpleName();
    }
}
