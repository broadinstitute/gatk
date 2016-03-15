package org.broadinstitute.hellbender.tools.spark.linkedreads;

import org.apache.avro.test.Simple;
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
import org.broadinstitute.hellbender.utils.read.GATKRead;
import scala.Tuple2;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

@CommandLineProgramProperties(
        summary = "Computes barcode overlap between genomic windows",
        oneLineSummary = "ComputeBarcodeOverlap",
        programGroup = SparkProgramGroup.class
)
public class ComputeBarcodeOverlapSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;

    @Argument(doc = "uri for the output file",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            optional = true)
    public String out;

    @Argument(doc = "Window size",
            shortName = "windowSize", fullName = "windowSize",
            optional = false)
    public Integer windowsSize = 10000;

    @Override
    public boolean requiresReads() {
        return true;
    }

    @Override
    public ReadFilter makeReadFilter() {
        return super.makeReadFilter()
                .and(read -> !read.isUnmapped())
                .and(read -> !read.failsVendorQualityCheck())
                .and(GATKRead::isFirstOfPair)
                .and(read -> !read.isDuplicate())
                .and(read -> !read.isSecondaryAlignment())
                .and(read -> !read.isSupplementaryAlignment())
                .and(read -> read.hasAttribute("BX"));
    }

    @Override
    protected void runTool(final JavaSparkContext ctx) {
        final JavaRDD<GATKRead> reads = getReads();


        final JavaPairRDD<String, Set<SimpleInterval>> windowsByBarcode = reads
                .mapToPair(read -> new Tuple2<>(read.getAttributeAsString("BX"), new SimpleInterval(read.getContig(), read.getStart() - read.getStart() % windowsSize + 1, read.getStart() + (windowsSize - read.getStart() % windowsSize))))
                .aggregateByKey(
                new HashSet<>(),
                (windowSet, window) -> {
                    windowSet.add(window);
                    return windowSet;
                },
                (barcodeSet1, barcodeSet2) -> {
                    barcodeSet1.addAll(barcodeSet2);
                    return barcodeSet1;
                });
        final JavaPairRDD<Tuple2<SimpleInterval, SimpleInterval>, Integer> pairCounts = windowsByBarcode.flatMapToPair(windowsForBarcode -> {
            final List<Tuple2<Tuple2<SimpleInterval, SimpleInterval>, Integer>> results = new ArrayList<>();
            final Set<SimpleInterval> windows = windowsForBarcode._2;
            for (final SimpleInterval window1 : windows) {
                for (final SimpleInterval window2 : windows) {
                    if (window1.getContig().compareTo(window2.getContig()) > 0) {
                        continue;
                    }
                    if (window1.getContig().equals(window2.getContig()) && window1.getStart() > window2.getStart()) {
                        continue;
                    }
                    results.add(new Tuple2<>(new Tuple2<>(window1, window2), 1));
                }
            }
            return results;
        }).reduceByKey((value1, value2) -> value1 + value2);

        pairCounts.saveAsTextFile(out);
    }
}
