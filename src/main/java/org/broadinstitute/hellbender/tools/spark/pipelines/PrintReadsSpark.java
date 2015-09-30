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
import org.broadinstitute.hellbender.cmdline.programgroups.SparkProgramGroup;
import org.broadinstitute.hellbender.engine.spark.SparkCommandLineProgram;
import org.broadinstitute.hellbender.engine.spark.datasources.ReadsSparkSink;
import org.broadinstitute.hellbender.engine.spark.datasources.ReadsSparkSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadsWriteFormat;

import java.io.IOException;
import java.util.List;

@CommandLineProgramProperties(summary = "Print reads from the input BAM", oneLineSummary = "Print reads from the input BAM", programGroup = SparkProgramGroup.class)
public final class PrintReadsSpark extends SparkCommandLineProgram {

    public static final String SHARDED_OUTPUT_LONG_ARG = "shardedOutput";
    public static final String SHARDED_OUTPUT_SHORT_ARG = "shardedOutput";

    private static final long serialVersionUID = 1L;

    @ArgumentCollection
    public ReadInputArgumentCollection readArguments= new RequiredReadInputArgumentCollection();

    @ArgumentCollection
    public IntervalArgumentCollection intervalArgumentCollection = new OptionalIntervalArgumentCollection();

    @Argument(doc = "uri for the output file: a local file path",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            optional = false)
    public String output;

    @Argument(doc = "If specified, shard the output bam", shortName = SHARDED_OUTPUT_SHORT_ARG, fullName = SHARDED_OUTPUT_LONG_ARG, optional = true)
    public boolean shardedOutput = false;

    @Override
    protected void runPipeline(final JavaSparkContext ctx) {
        if ( readArguments.getReadFilesNames().size() != 1 ) {
            throw new UserException("This tool only accepts a single bam/sam/cram as input");
        }
        final String bam = readArguments.getReadFilesNames().get(0);

        final ReadsSparkSource readSource = new ReadsSparkSource(ctx);
        final SAMFileHeader readsHeader = ReadsSparkSource.getHeader(ctx, bam);

        if (readsHeader.getSortOrder() != SAMFileHeader.SortOrder.coordinate){
            //https://github.com/broadinstitute/hellbender/issues/929
            throw new UserException("PrintReadsSpark: Only coordinate-sorted files are currently supported");
        }

        /*
         * If no intervals are given, we want all reads, mapped and unmapped.
         */
        final List<SimpleInterval> intervals =  intervalArgumentCollection.intervalsSpecified() ?
                intervalArgumentCollection.getIntervals(readsHeader.getSequenceDictionary()) :
                null;
        final JavaRDD<GATKRead> reads= readSource.getParallelReads(bam, intervals);

        try {
            ReadsSparkSink.writeReads(ctx, output, reads, readsHeader, shardedOutput ? ReadsWriteFormat.SHARDED : ReadsWriteFormat.SINGLE);
        } catch (final IOException e) {
            throw new GATKException("unable to write bam: " + e);
        }
    }
}