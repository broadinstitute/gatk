package org.broadinstitute.hellbender.tools.spark.linkedreads;

import htsjdk.samtools.util.IntervalTree;
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

    @Argument(doc = "uri for the input linked reads file",
            shortName = "linked reads file", fullName = "Linked reads file created by CollectLinkedReadCoverage tool",
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

        final JavaPairRDD<SimpleInterval, Iterable<String>> intervalsToLinkedReads = barcodeIntervals.join(barcodeLinkedReads).values().mapToPair(value -> value).groupByKey();

        intervalsToLinkedReads.saveAsTextFile(out);
    }

    private static List<SimpleInterval> getIntervalsForRead(final IntervalsSkipList<SimpleInterval> intervalsSkipList, final GATKRead read) {
        return intervalsSkipList.getOverlapping(new SimpleInterval(read));
    }
}
