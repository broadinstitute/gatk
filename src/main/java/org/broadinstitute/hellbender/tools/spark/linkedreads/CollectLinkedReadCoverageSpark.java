package org.broadinstitute.hellbender.tools.spark.linkedreads;

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
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import scala.Tuple2;

import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
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

    @Override
    public boolean requiresReads() { return true; }

    @Override
    public boolean requiresReference() { return true; }

    @Override
    protected void runTool(final JavaSparkContext ctx) {
        final JavaRDD<GATKRead> reads = getReads();

        final JavaPairRDD<String, Map<String, IntervalTree<Integer>>> intervalsByBarcode =
                reads.mapToPair(read -> new Tuple2<>(read.getAttributeAsString("BX"), read))
                        .aggregateByKey(
                                new HashMap<>(),
                                CollectLinkedReadCoverageSpark::addReadToIntervals,
                                CollectLinkedReadCoverageSpark::combineIntervalLists
                        );

        intervalsByBarcode.saveAsTextFile(out);
    }

    @VisibleForTesting
    static Map<String, IntervalTree<Integer>> addReadToIntervals(Map<String, IntervalTree<Integer>> intervalList, GATKRead read) {
        final Interval sloppedReadInterval = new Interval(read.getContig(), read.getStart() - 1000, read.getStart()+1000);
        if (! intervalList.containsKey(read.getContig())) {
            intervalList.put(read.getContig(), new IntervalTree<>());
        }
        final IntervalTree<Integer> contigIntervals = intervalList.get(sloppedReadInterval.getContig());

        final Iterator<IntervalTree.Node<Integer>> iterator = contigIntervals.overlappers(sloppedReadInterval.getStart(), sloppedReadInterval.getEnd());
        int start = read.getStart();
        int end = read.getEnd();
        int value = 1;
        if (iterator.hasNext()) {
            final IntervalTree.Node<Integer> next = iterator.next();
            final int currentStart = next.getStart();
            final int currentEnd = next.getStart();
            final int currentValue = next.getValue();
            start = Math.min(currentStart, read.getStart());
            end = Math.max(currentEnd, read.getEnd());
            value = currentValue + 1;
            iterator.remove();
        }
        while (iterator.hasNext()) {
            final IntervalTree.Node<Integer> next = iterator.next();
            final int currentEnd = next.getStart();
            final int currentValue = next.getValue();
            end = currentEnd;
            value = value + currentValue;
            iterator.remove();
        }
        contigIntervals.put(start, end, value);
        return intervalList;
    }

    static Map<String, IntervalTree<Integer>> combineIntervalLists(Map<String, IntervalTree<Integer>> intervalList1, Map<String, IntervalTree<Integer>> intervalList2) {
        Map<String, IntervalTree<Integer>> combinedList = new HashMap<>();
        Stream.concat(intervalList1.keySet().stream(), intervalList2.keySet().stream()).
                forEach(contigName -> combinedList.put(contigName, mergeIntervalTrees(intervalList1.get(contigName), intervalList2.get(contigName))));
        return intervalList1;

    }

    @VisibleForTesting
    static IntervalTree<Integer> mergeIntervalTrees(final IntervalTree<Integer> tree1, final IntervalTree<Integer> tree2) {

        if (tree1 == null || tree1.size() == 0) return tree2;
        if (tree2 == null || tree2.size() == 0) return tree1;

        tree1.iterator().forEachRemaining(node -> tree2.put(node.getStart(), node.getEnd(), node.getValue()));

        final IntervalTree<Integer> mergedTree = new IntervalTree<>();

        while (tree2.size() > 0) {
            final IntervalTree.Node<Integer> current = tree2.min();
            final int currentStart = current.getStart();
            int currentEnd = current.getEnd();
            int currentValue = current.getValue();

            final Iterator<IntervalTree.Node<Integer>> overlappers = tree2.overlappers(current.getStart(), current.getEnd());
            while (overlappers.hasNext()) {
                IntervalTree.Node<Integer> overlapper = overlappers.next();
                if (overlapper == current) {
                    continue;
                }
                currentEnd = overlapper.getEnd();
                currentValue = currentValue + overlapper.getValue();
                overlappers.remove();
            }
            mergedTree.put(currentStart, currentEnd, currentValue);
            tree2.remove(current.getStart(), current.getEnd());
        }
        return mergedTree;
    }
}
