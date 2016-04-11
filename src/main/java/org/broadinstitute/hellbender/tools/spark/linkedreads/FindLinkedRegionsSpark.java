package org.broadinstitute.hellbender.tools.spark.linkedreads;

import htsjdk.samtools.util.IntervalTree;
import org.apache.spark.HashPartitioner;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.SparkProgramGroup;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.collections.IntervalsSkipList;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import scala.Tuple2;

import java.util.*;
import java.util.stream.Collectors;

@CommandLineProgramProperties(
        summary = "Computes barcode overlap between genomic windows",
        oneLineSummary = "ComputeBarcodeOverlap",
        programGroup = SparkProgramGroup.class
)
public class FindLinkedRegionsSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;
    public static final int COVERAGE_BIN_WIDTH = 1000;

    @Argument(doc = "Linked reads file created by CollectLinkedReadCoverage tool",
            shortName = "linkedReadsFile", fullName = "linkedReadsFile",
            optional = false)
    public String linkedReadsFile;

    @Argument(doc = "uri for the output file",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            optional = true)
    public String out;

    @Override
    public boolean requiresReads() {
        return true;
    }

    @Override
    public boolean requiresIntervals() {
        return true;
    }

    @Override
    public ReadFilter makeReadFilter() {
        return super.makeReadFilter().and(new LinkedReadAnalysisFilter());
    }

    @Override
    protected void runTool(final JavaSparkContext ctx) {
        final List<SimpleInterval> intervalList = getIntervals();
        final IntervalsSkipList<SimpleInterval> intervalsSkipList = new IntervalsSkipList<>(intervalList);

        final JavaRDD<GATKRead> reads = getReads();

        final JavaPairRDD<String, SimpleInterval> barcodeIntervals = reads.flatMapToPair(read -> {
            final List<SimpleInterval> overlappingIntervals = getIntervalsForRead(intervalsSkipList, read);
            return overlappingIntervals
                    .stream()
                    .map(interval -> new Tuple2<>(read.getAttributeAsString("BX"), interval))
                    .collect(Collectors.toList());
        }).distinct();


        final JavaRDD<String> linkedReadsBedFile = ctx.textFile(linkedReadsFile);

        final JavaPairRDD<String, String> barcodeLinkedReads = linkedReadsBedFile.mapToPair(line -> new Tuple2<>(line.split("\t")[3], line));

        final JavaPairRDD<String, Tuple2<SimpleInterval, String>> intervalsAndLinkedReadJoinedByBarcodes = barcodeIntervals.join(barcodeLinkedReads);

        final JavaPairRDD<Tuple2<SimpleInterval, SimpleInterval>, Integer> coverageCountsByIntervalAndBin = intervalsAndLinkedReadJoinedByBarcodes.flatMapToPair(join -> {
            final Tuple2<SimpleInterval, String> intervalLinkedReadPair = join._2();
            final List<Tuple2<Tuple2<SimpleInterval, SimpleInterval>, Integer>> results = new ArrayList<>();
            final SimpleInterval linkedReadInterval = getLinkedReadInterval(intervalLinkedReadPair._2());
            for (int binStart = linkedReadInterval.getStart() - linkedReadInterval.getStart() % COVERAGE_BIN_WIDTH; binStart < linkedReadInterval.getEnd() + COVERAGE_BIN_WIDTH; binStart += COVERAGE_BIN_WIDTH) {
                final SimpleInterval bin = new SimpleInterval(linkedReadInterval.getContig(), binStart, binStart + COVERAGE_BIN_WIDTH - 1);
                results.add(new Tuple2<>(new Tuple2<>(intervalLinkedReadPair._1(), bin), 1));
            }
            return results;
        }).reduceByKey((count1, count2) -> count1 + count2);

        final JavaPairRDD<SimpleInterval, Tuple2<SimpleInterval, Integer>> intervalBinCoverages =
                coverageCountsByIntervalAndBin.mapToPair(intervalBinCovCount -> new Tuple2<>(intervalBinCovCount._1._1, new Tuple2<>(intervalBinCovCount._1._2, intervalBinCovCount._2)))
                .partitionBy(new HashPartitioner(intervalList.size()));

        intervalBinCoverages.saveAsTextFile(out);
    }

    private SimpleInterval getLinkedReadInterval(final String linkedReadString) {
        final String[] fields = linkedReadString.split("\t");
        return new SimpleInterval(fields[0], Integer.valueOf(fields[1]), Integer.valueOf(fields[2]));
    }

    private static List<SimpleInterval> getIntervalsForRead(final IntervalsSkipList<SimpleInterval> intervalsSkipList, final GATKRead read) {
        return intervalsSkipList.getOverlapping(new SimpleInterval(read));
    }
}
