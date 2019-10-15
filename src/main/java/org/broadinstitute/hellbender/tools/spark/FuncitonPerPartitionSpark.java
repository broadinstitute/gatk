package org.broadinstitute.hellbender.tools.spark;

import htsjdk.samtools.SAMFileHeader;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.spark.SparkUtils;
import picard.cmdline.programgroups.ReadDataManipulationProgramGroup;

import java.util.Collections;
import java.util.List;

/**
 * A quick dummy class for doing per-readgroup operations in spark. If you need to mutate the reads based on the results
 * refer to MarkDuplicatesSpark (or ask) and the operation can be setup (as its inefficient if not done right).
 */
@DocumentedFeature
@CommandLineProgramProperties(
        summary ="Marks duplicates on spark",
        oneLineSummary ="MarkDuplicates on Spark",
        programGroup = ReadDataManipulationProgramGroup.class)
public final class MarkDuplicatesSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;

    @Override
    public boolean requiresReads() { return true; }

    @Argument(doc = "the output bam", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, optional = false)
    protected String output;

    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        return Collections.singletonList(ReadFilterLibrary.ALLOW_ALL_READS);
    }

    @Override
    public ReadInputMergingPolicy getReadInputMergingPolicy() {
        return ReadInputMergingPolicy.concatMerge;
    }

    /**
     * MAIN METHOD TO FILL IN YOUR READ BEHAVIOR
     *
     * this method should return an RDD of whatever your return type is per each readgroup.
     */
    public static JavaRDD<String> processReads(final JavaPairRDD<String, Iterable<GATKRead>> readsByName, final SAMFileHeader header) {
        // inside this mapping operation is a TUPLE2<String (readname), Iterable<GATKRead> (reads with that name)>
        return readsByName.map(pair -> pair._1);
    }

    @Override
    protected void runTool(final JavaSparkContext ctx) {
        final SAMFileHeader mergedHeader = getHeaderForReads();

        JavaRDD<GATKRead> reads = getReads();
        // If the input isn't queryname sorted, sort it before duplicate marking
        final JavaRDD<GATKRead> sortedReadsForWork = SparkUtils.querynameSortReadsIfNecessary(reads, numReducers, getHeaderForReads());

        JavaRDD<String> processedReads = processReads(RevertSamSpark.spanReadsByKey(sortedReadsForWork), mergedHeader);

        //FILL IN YOUR OWN RETURN RESULTS
    }

}
