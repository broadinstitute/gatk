package org.broadinstitute.hellbender.tools.spark.pipelines;

import htsjdk.samtools.SAMFileHeader;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.utils.spark.SparkUtils;
import picard.cmdline.programgroups.ReadDataManipulationProgramGroup;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.Collections;
import java.util.List;

@DocumentedFeature
@CommandLineProgramProperties(summary = "Sorts the input SAM/BAM/CRAM",
        oneLineSummary = "SortSam on Spark (works on SAM/BAM/CRAM)",
        programGroup = ReadDataManipulationProgramGroup.class)
@BetaFeature
public final class SortSamSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;

    public static final String SORT_ORDER_LONG_NAME = "sort-order";

    @Override
    public boolean requiresReads() { return true; }

    @Argument(doc="the output file path", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, optional = false)
    private String outputFile;

    @Argument(doc="sort order of the output file", fullName = SORT_ORDER_LONG_NAME, optional = true)
    private SparkSortOrder sortOrder = SparkSortOrder.coordinate;

    /**
     * SortOrders that have corresponding implementations for spark.
     * These correspond to a subset of {@link SAMFileHeader.SortOrder}.
     */
    private enum SparkSortOrder {
        coordinate(SAMFileHeader.SortOrder.coordinate),
        queryname(SAMFileHeader.SortOrder.queryname);

        private final SAMFileHeader.SortOrder order;

        SparkSortOrder(SAMFileHeader.SortOrder order) {
            this.order = order;
        }

        public SAMFileHeader.SortOrder getSamOrder() {
             return order;
        }
    }

    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        return Collections.singletonList(ReadFilterLibrary.ALLOW_ALL_READS);
    }

    @Override
    protected void onStartup() {
        super.onStartup();
    }

    @Override
    protected void runTool(final JavaSparkContext ctx) {
        final JavaRDD<GATKRead> reads = getReads();
        final int numReducers = getRecommendedNumReducers();
        logger.info("Using %d reducers", numReducers);

        final SAMFileHeader header = getHeaderForReads();
        header.setSortOrder(sortOrder.getSamOrder());

        final JavaRDD<GATKRead> readsToWrite;
        if (shardedOutput) {
            readsToWrite = SparkUtils.sortReadsAccordingToHeader(reads, header, numReducers);
        } else {
            readsToWrite = reads;
        }
        writeReads(ctx, outputFile, readsToWrite, header);
    }
}
