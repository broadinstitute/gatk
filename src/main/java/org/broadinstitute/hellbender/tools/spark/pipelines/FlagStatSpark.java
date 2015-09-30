package org.broadinstitute.hellbender.tools.spark.pipelines;

import htsjdk.samtools.SAMFileHeader;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.ArgumentCollection;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.IntervalArgumentCollection;
import org.broadinstitute.hellbender.cmdline.argumentcollections.OptionalIntervalArgumentCollection;
import org.broadinstitute.hellbender.cmdline.argumentcollections.ReadInputArgumentCollection;
import org.broadinstitute.hellbender.cmdline.argumentcollections.RequiredReadInputArgumentCollection;
import org.broadinstitute.hellbender.cmdline.programgroups.DataFlowProgramGroup;
import org.broadinstitute.hellbender.cmdline.programgroups.SparkProgramGroup;
import org.broadinstitute.hellbender.engine.spark.SparkCommandLineProgram;
import org.broadinstitute.hellbender.engine.spark.datasources.ReadsSparkSource;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.FlagStat;
import org.broadinstitute.hellbender.tools.FlagStat.FlagStatus;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.dataflow.BucketUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.List;

@CommandLineProgramProperties(summary ="runs FlagStat on Spark", oneLineSummary = "FlagStat", programGroup = SparkProgramGroup.class)
public final class FlagStatSpark extends SparkCommandLineProgram {

    private static final long serialVersionUID = 1L;

    @ArgumentCollection
    public ReadInputArgumentCollection readArguments= new RequiredReadInputArgumentCollection();;

    @ArgumentCollection
    public IntervalArgumentCollection intervalArgumentCollection = new OptionalIntervalArgumentCollection();

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
        SAMFileHeader readsHeader = ReadsSparkSource.getHeader(ctx, bam);

        /**
         * If no intervals specified - use all reads, mapped and unmapped.
         */
        final List<SimpleInterval> intervals = intervalArgumentCollection.intervalsSpecified() ? intervalArgumentCollection.getIntervals(readsHeader.getSequenceDictionary())
                : null;

        final JavaRDD<GATKRead> reads = readSource.getParallelReads(bam, intervals);
        final FlagStatus result = reads.aggregate(new FlagStatus(), FlagStatus::add, FlagStatus::merge);
        final File file = new File(out);
        System.out.println(result);
        try(final OutputStream outputStream = BucketUtils.createFile(file.getPath(), null);
            final PrintStream ps = new PrintStream(outputStream)) {
            ps.print(result);
        } catch(final IOException e){
            throw new UserException.CouldNotCreateOutputFile(file, e);
        }
    }
}
