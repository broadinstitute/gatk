package org.broadinstitute.hellbender.tools.spark.pipelines;

import htsjdk.samtools.SAMFileHeader;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.SparkProgramGroup;
import org.broadinstitute.hellbender.engine.spark.SparkCommandLineProgram;
import org.broadinstitute.hellbender.engine.spark.datasources.ReadsSparkSink;
import org.broadinstitute.hellbender.engine.spark.datasources.ReadsSparkSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadCoordinateComparator;
import org.broadinstitute.hellbender.utils.read.ReadsWriteFormat;
import scala.Tuple2;

import java.io.IOException;

@CommandLineProgramProperties(summary = "Sorts the input BAM", oneLineSummary = "Sorts a BAM file", programGroup = SparkProgramGroup.class)
public final class SortBamSpark extends SparkCommandLineProgram {

    private static final long serialVersionUID = 1L;

    @Argument(doc="the output file path", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, optional = false)
    protected String outputFile;

    @Argument(doc = "uri for the input bam, either a local file path, hdfs:// path, or a gs:// bucket path",
            shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME, fullName = StandardArgumentDefinitions.INPUT_LONG_NAME,
            optional = false)
    protected String bam;

    @Argument(doc="the output parallelism, sets the number of reducers", shortName = "P", fullName = "parallelism", optional = true)
    protected int parallelism = 0;

    @Override
    protected void runPipeline(final JavaSparkContext ctx) {
        ReadsSparkSource readSource = new ReadsSparkSource(ctx);
        SAMFileHeader readsHeader = ReadsSparkSource.getHeader(ctx, bam);
        JavaRDD<GATKRead> reads = readSource.getParallelReads(bam, null);
        if (parallelism == 0) {
            parallelism = reads.partitions().size();
        }
        System.out.println("Using parallelism of " + parallelism);
        ReadCoordinateComparator comparator = new ReadCoordinateComparator(readsHeader);
        JavaRDD<GATKRead> sortedReads = reads
                .mapToPair(read -> new Tuple2<>(read, null))
                .sortByKey(comparator, true, parallelism)
                .keys();
        try {
            readsHeader.setSortOrder(SAMFileHeader.SortOrder.coordinate);
            ReadsSparkSink.writeReads(
                    ctx, outputFile, sortedReads, readsHeader, parallelism == 1 ? ReadsWriteFormat.SINGLE : ReadsWriteFormat.SHARDED);
        } catch (IOException e) {
            throw new GATKException("unable to write bam: " + e);
        }
    }

    @Override
    protected String getProgramName() {
        return getClass().getSimpleName();
    }
}
