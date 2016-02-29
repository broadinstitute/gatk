package org.broadinstitute.hellbender.tools.spark.linkedreads;

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

import java.util.HashSet;
import java.util.Set;

@CommandLineProgramProperties(
        summary = "Computes the coverage by both short and long reads in bins across the genome",
        oneLineSummary = "CollectBarcodedCoverage on Spark",
        programGroup = SparkProgramGroup.class
)
public class CollectBarcodeCoverageCountsSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;

    @Argument(doc = "uri for the output file",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            optional = false)
    public String out;

    @Argument(doc = "bin size",
            shortName = "binSize", fullName = "binSize",
            optional = false)
    public Integer binSize;

    @Override
    public boolean requiresReads() { return true; }

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

        final JavaPairRDD<SimpleInterval, Tuple2<Integer, Integer>> intervalCounts =
                reads.mapToPair(read -> new Tuple2<>(new SimpleInterval(read.getContig(), read.getStart() - read.getStart() % binSize + 1, read.getStart() + (binSize - read.getStart() % binSize)), read))
                        .aggregateByKey(
                                new Tuple2<Integer, Set<String>>(1, new HashSet<>()),
                                (aggregator, read) -> {
                                    int newCount = aggregator._1 + 1;
                                    Set<String> newBarcodeSet = aggregator._2;
                                    newBarcodeSet.add(read.getAttributeAsString("BX"));
                                    return new Tuple2<>(newCount, newBarcodeSet);
                                },
                                (aggregator1, aggregator2) -> {
                                    int newCount = aggregator1._1 + aggregator2._1;
                                    Set<String> newBarcodeSet = aggregator1._2;
                                    newBarcodeSet.addAll(aggregator2._2);
                                    return new Tuple2<>(newCount, newBarcodeSet);
                                }
                        )
                .mapValues(aggregator -> new Tuple2<>(aggregator._1, aggregator._2.size()));

        intervalCounts.sortByKey().coalesce(1).saveAsTextFile(out);
    }
}
