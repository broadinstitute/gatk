package org.broadinstitute.hellbender.tools.spark.linkedreads;

import com.google.api.services.genomics.model.Read;
import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.util.Interval;
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
import org.broadinstitute.hellbender.utils.read.GATKRead;
import scala.Tuple2;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

@CommandLineProgramProperties(
        summary = "Computes the coverage by long molecules from linked-read data",
        oneLineSummary = "CollectLinkedReadCoverage on Spark",
        programGroup = SparkProgramGroup.class
)
public class CollectLinkedReadCoverageSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;

    @Argument(doc = "uri for the output file",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            optional = true)
    public String out;

    @Argument(doc = "cluster size",
            shortName = "clusterSize", fullName = "clusterSize",
            optional = true)
    public int clusterSize = 5000;

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

        final JavaPairRDD<String, Map<String, IntervalTree<List<ReadInfo>>>> barcodeIntervals =
                reads.mapToPair(read -> new Tuple2<>(read.getAttributeAsString("BX"), read))
                .aggregateByKey(
                        new HashMap<>(),
                        (aggregator, read) -> addReadToIntervals(aggregator, read, clusterSize),
                        (intervalTree1, intervalTree2) -> combineIntervalLists(intervalTree1, intervalTree2, clusterSize)
                );
        final JavaRDD<String> intervalsByBarcode =
                barcodeIntervals.flatMap(CollectLinkedReadCoverageSpark::intervalTreeToString);


        intervalsByBarcode.saveAsTextFile(out);
    }

    static List<String> intervalTreeToString(final Tuple2<String, Map<String, IntervalTree<List<ReadInfo>>>> barcodeLocations) {
        final List<String> outputLines = new ArrayList<>();
        final String barcode = barcodeLocations._1;
        final Map<String, IntervalTree<List<ReadInfo>>> stringIntervalTreeMap = barcodeLocations._2;
        for (final String contig : stringIntervalTreeMap.keySet()) {
            final IntervalTree<List<ReadInfo>> contigIntervalTree = stringIntervalTreeMap.get(contig);
            final Iterator<IntervalTree.Node<List<ReadInfo>>> nodeIterator = contigIntervalTree.iterator();
            while (nodeIterator.hasNext()) {
                final StringBuilder out = new StringBuilder();
                final IntervalTree.Node<List<ReadInfo>> node = nodeIterator.next();
                out.append(contig);
                out.append("\t");
                out.append(node.getStart());
                out.append("\t");
                out.append(node.getEnd());
                out.append("\t");
                out.append(barcode);
                out.append("\t");
                out.append(node.getValue().size());
                out.append("\t");
                out.append("+");
                out.append("\t");
                out.append(node.getStart());
                out.append("\t");
                out.append(node.getEnd());
                out.append("\t");
                out.append("0,0,255");
                out.append("\t");
                out.append(node.getValue().size());
                List<ReadInfo> results = node.getValue();
                results.sort((o1, o2) -> new Integer(o1.start).compareTo(o2.start));
                out.append("\t");
                out.append(results.stream().map(r -> String.valueOf(r.end - r.start)).collect(Collectors.joining(",")));
                out.append("\t");
                out.append(results.stream().map(r -> String.valueOf(r.start - node.getStart())).collect(Collectors.joining(",")));
                outputLines.add(out.toString());
            }
        }

        return outputLines;
    }

    @VisibleForTesting
    static Map<String, IntervalTree<List<ReadInfo>>> addReadToIntervals(final Map<String, IntervalTree<List<ReadInfo>>> intervalList, final GATKRead read, final int clusterSize) {
        final Interval sloppedReadInterval = new Interval(read.getContig(), read.getStart() - clusterSize, read.getStart()+clusterSize);
        if (! intervalList.containsKey(read.getContig())) {
            intervalList.put(read.getContig(), new IntervalTree<>());
        }
        final IntervalTree<List<ReadInfo>> contigIntervals = intervalList.get(sloppedReadInterval.getContig());

        final Iterator<IntervalTree.Node<List<ReadInfo>>> iterator = contigIntervals.overlappers(sloppedReadInterval.getStart(), sloppedReadInterval.getEnd());
        int start = read.getStart();
        int end = read.getEnd();
        List<ReadInfo> value = new ArrayList<>();
        final ReadInfo readInfo = new ReadInfo(read.getStart(), read.getEnd());
        value.add(readInfo);
        if (iterator.hasNext()) {
            final IntervalTree.Node<List<ReadInfo>> next = iterator.next();
            final int currentStart = next.getStart();
            final int currentEnd = next.getStart();
            final List<ReadInfo> currentValue = next.getValue();
            start = Math.min(currentStart, read.getStart());
            end = Math.max(currentEnd, read.getEnd());
            value.addAll(currentValue);
            iterator.remove();
        }
        while (iterator.hasNext()) {
            final IntervalTree.Node<List<ReadInfo>> next = iterator.next();
            final int currentEnd = next.getStart();
            final List<ReadInfo> currentValue = next.getValue();
            end = currentEnd;
            value.addAll(currentValue);
            iterator.remove();
        }
        contigIntervals.put(start, end, value);
        return intervalList;
    }

    static Map<String, IntervalTree<List<ReadInfo>>> combineIntervalLists(final Map<String, IntervalTree<List<ReadInfo>>> intervalList1,
                                                                   final Map<String, IntervalTree<List<ReadInfo>>> intervalList2,
                                                                   final int clusterSize) {
        final Map<String, IntervalTree<List<ReadInfo>>> combinedList = new HashMap<>();
        Stream.concat(intervalList1.keySet().stream(), intervalList2.keySet().stream()).
                forEach(contigName -> combinedList.put(contigName, mergeIntervalTrees(intervalList1.get(contigName), intervalList2.get(contigName), clusterSize)));
        return combinedList;

    }

    @VisibleForTesting
    static IntervalTree<List<ReadInfo>> mergeIntervalTrees(final IntervalTree<List<ReadInfo>> tree1, final IntervalTree<List<ReadInfo>> tree2, final int clusterSize) {

        if (tree1 == null || tree1.size() == 0) return tree2;
        if (tree2 == null || tree2.size() == 0) return tree1;

        tree1.iterator().forEachRemaining(node -> tree2.put(node.getStart(), node.getEnd(), node.getValue()));

        final IntervalTree<List<ReadInfo>> mergedTree = new IntervalTree<>();

        while (tree2.size() > 0) {
            final IntervalTree.Node<List<ReadInfo>> current = tree2.min();
            final int currentStart = current.getStart();
            int currentEnd = current.getEnd();
            List<ReadInfo> currentValue = current.getValue();

            final Iterator<IntervalTree.Node<List<ReadInfo>>> overlappers = tree2.overlappers(current.getStart() - clusterSize, current.getEnd() + clusterSize);
            while (overlappers.hasNext()) {
                final IntervalTree.Node<List<ReadInfo>> overlapper = overlappers.next();
                if (overlapper == current) {
                    continue;
                }
                currentEnd = overlapper.getEnd();
                currentValue.addAll(overlapper.getValue());
                overlappers.remove();
            }
            mergedTree.put(currentStart, currentEnd, currentValue);
            tree2.remove(current.getStart(), current.getEnd());
        }
        return mergedTree;
    }

    static class ReadInfo {
        public ReadInfo(final int start, final int end) {
            this.start = start;
            this.end = end;
        }

        int start;
        int end;
    }
}
