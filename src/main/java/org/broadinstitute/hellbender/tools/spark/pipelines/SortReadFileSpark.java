package org.broadinstitute.hellbender.tools.spark.pipelines;

import htsjdk.samtools.SAMFileHeader;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.SparkProgramGroup;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadCoordinateComparator;
import scala.Tuple2;

@CommandLineProgramProperties(summary = "Sorts the input SAM/BAM/CRAM",
        oneLineSummary = "SortSam on Spark (works on SAM/BAM/CRAM)",
        programGroup = SparkProgramGroup.class)
public final class SortReadFileSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;

    @Override
    public boolean requiresReads() { return true; }

    @Argument(doc="the output file path", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, optional = false)
    protected String outputFile;

    @Override
    public ReadFilter makeReadFilter() {
        return ReadFilterLibrary.ALLOW_ALL_READS;
    }

    @Override
    protected void runTool(final JavaSparkContext ctx) {
        JavaRDD<GATKRead> reads = getReads();
        int numReducers = getRecommendedNumReducers();
        logger.info("Using %s reducers" + numReducers);

        final SAMFileHeader readsHeader = getHeaderForReads();
        ReadCoordinateComparator comparator = new ReadCoordinateComparator(readsHeader);
        JavaRDD<GATKRead> sortedReads;
        if (shardedOutput) {
            sortedReads = reads
                    .mapToPair(read -> new Tuple2<>(read, null))
                    .sortByKey(comparator, true, numReducers)
                    .keys();
        } else {
            sortedReads = reads; // sorting is done by writeReads below
        }
        readsHeader.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        writeReads(ctx, outputFile, sortedReads);
    }
}
