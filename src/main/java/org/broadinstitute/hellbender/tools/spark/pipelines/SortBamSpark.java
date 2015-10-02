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
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadCoordinateComparator;
import org.broadinstitute.hellbender.utils.read.ReadsWriteFormat;
import scala.Tuple2;

import java.io.IOException;

@CommandLineProgramProperties(summary = "Sorts the input BAM", oneLineSummary = "Sorts a BAM file", programGroup = SparkProgramGroup.class)
public final class SortBamSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;

    @Override
    public boolean requiresReads() { return true; }

    @Argument(doc="the output file path", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, optional = false)
    protected String outputFile;

    @Argument(doc="The output parallelism, sets the number of reducers. Defaults to the number of partitions in the input.",
            shortName = "P", fullName = "parallelism", optional = true)
    protected int parallelism = 0;

    @Override
    protected void runTool(final JavaSparkContext ctx) {
        JavaRDD<GATKRead> reads = getReads();

        if (parallelism == 0) { // use the number of partitions in the input
            parallelism = reads.partitions().size();
        }
        System.out.println("Using parallelism of " + parallelism);

        final SAMFileHeader readsHeader = getHeaderForReads();
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
}
