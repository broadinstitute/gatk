package org.broadinstitute.hellbender.tools.spark.linkedreads;

import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.SparkProgramGroup;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import scala.Tuple2;

import java.util.HashSet;
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
    protected void runTool(final JavaSparkContext ctx) {
        final JavaRDD<GATKRead> reads = getReads();


        final JavaPairRDD<SimpleInterval, Set<String>> barcodesByWindow = reads.mapToPair(read -> new Tuple2<>(new SimpleInterval(read.getContig(), read.getStart() - read.getStart() % windowsSize + 1, read.getStart() + (windowsSize - read.getStart() % windowsSize)), read)).aggregateByKey(
                new HashSet<>(),
                (barcodeSet, read) -> {
                    final Set<String> newBarcodeSet = new HashSet<>(barcodeSet);
                    newBarcodeSet.add(read.getAttributeAsString("BX"));
                    return newBarcodeSet;
                },
                (barcodeSet1, barcodeSet2) -> {
                    final Set<String> newBarcodeSet = new HashSet<>(barcodeSet1);
                    newBarcodeSet.addAll(barcodeSet2);
                    return newBarcodeSet;
                });
        final JavaPairRDD<Tuple2<SimpleInterval, SimpleInterval>, Integer> barcodeOverlapsByWindowPairs = barcodesByWindow.cartesian(barcodesByWindow).filter(windowPair -> {
            // only keep windows where the first window comes before the second
            final SimpleInterval interval1 = windowPair._1._1;
            final SimpleInterval interval2 = windowPair._2._1;
            if (interval1.getContig().compareTo(interval2.getContig()) > 0) {
                return false;
            }
            if (interval1.getContig().compareTo(interval2.getContig()) < 0) {
                return true;
            }
            if (interval1.getStart() < interval2.getStart()) {
                return true;
            }
            return false;
        }).mapToPair(windowPair -> {
            final Tuple2<SimpleInterval, SimpleInterval> newkey = new Tuple2<>(windowPair._1._1, windowPair._2._1);
            final Set<String> combinedSet = new HashSet<>(windowPair._1()._2());
            combinedSet.retainAll(windowPair._2._2);
            return new Tuple2<>(newkey, combinedSet.size());
        });


        barcodeOverlapsByWindowPairs.saveAsObjectFile(out);
    }
}
