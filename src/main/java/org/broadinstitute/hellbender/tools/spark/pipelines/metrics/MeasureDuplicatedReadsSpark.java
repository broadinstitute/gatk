package org.broadinstitute.hellbender.tools.spark.pipelines.metrics;

import com.google.common.collect.Iterators;
import com.google.common.collect.PeekingIterator;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.api.java.function.FlatMapFunction2;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CoverageAnalysisProgramGroup;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/**
 * Calculate the overall number of reads in a SAM/BAM file
 *
 * <h3>Input</h3>
 * <ul>
 *     <li>A single BAM file</li>
 * </ul>
 *
 * <h3>Output</h3>
 * <ul>
 *     <li>A text file containing number of reads</li>
 * </ul>
 *
 * <h3>Example</h3>
 *
 * <h4>Output number of reads to file</h4>
 * <pre>
 *   gatk CountReadsSpark \
 *     -I input_reads.bam \
 *     -O read_count.txt
 * </pre>
 *
 * <h4>Print read count</h4>
 * <pre>
 *   gatk CountReadsSpark \
 *     -I input_reads.bam
 * </pre>
 */
@DocumentedFeature
@BetaFeature
@CommandLineProgramProperties(
        summary = "Counts reads in the input SAM/BAM",
        oneLineSummary = "Counts reads in the input SAM/BAM",
        programGroup = CoverageAnalysisProgramGroup.class
)
public final class MeasureDuplicatedReadsSpark extends GATKSparkTool {

    private static final long serialVersionUID = 1L;

    @Override
    public boolean requiresReads() { return true; }

    @Argument(
            doc = "uri for the output file: a local file path",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            optional = true
    )
    public String out;

    @Override
    protected void runTool(final JavaSparkContext ctx) {
        final JavaRDD<GATKRead> reads = getReads();

        final long count = reads.count();
        System.out.println(count);

        final int numParititions = reads.partitions().size();
        System.out.println("There are " + numParititions +" partitions of reads generated for this bam");
        final String firstNameInBam = reads.first().getName();
        // Find the first group in each partition
        List<List<GATKRead>> firstReadNameGroupInEachPartition = reads
                .mapPartitions(it -> { PeekingIterator<GATKRead> current = Iterators.peekingIterator(it);
                    List<GATKRead> firstGroup = new ArrayList<>(2);
                    boolean pastFirstGroup = false;
                    int i = 0;
                    firstGroup.add(current.next());
                    GATKRead last = firstGroup.get(0);
                    while (current.hasNext() && current.peek().getName().equals(last.getName())) {
                        if (last.equals(current.peek())) {
                            throw new GATKException(String.format("Read %s was duplicated in a partition at position &d in the partition", last, i));
                        }

                        if (!pastFirstGroup && current.peek().equals(last)) {
                            firstGroup.add(current.peek());
                        } else {
                            pastFirstGroup = true;
                        }
                        last = current.next();
                        i++;
                    }
                    return Iterators.singletonIterator(firstGroup);
                })
                .collect();

        // Checking for pathological cases (read name groups that span more than 2 partitions)
        String groupName = null;
        for (List<GATKRead> group : firstReadNameGroupInEachPartition) {
            if (group!=null && !group.isEmpty()) {
                // If a read spans multiple partitions we expect its name to show up multiple times and we don't expect this to work properly
                if (groupName != null && group.get(0).getName().equals(groupName)) {
                    throw new GATKException(String.format("The read name '%s' appeared across multiple partitions this could indicate there was a problem " +
                            "with the sorting or that the rdd has too many partitions, check that the file is queryname sorted and consider decreasing the number of partitions", groupName));
                }
                groupName =  group.get(0).getName();
            }
        }

        // Shift left, so that each partition will be joined with the first read group from the _next_ partition
        List<List<GATKRead>> firstReadInNextPartition = new ArrayList<>(firstReadNameGroupInEachPartition.subList(1, numParititions));
        firstReadInNextPartition.add(null); // the last partition does not have any reads to add to it

        // Join the reads with the first read from the _next_ partition, then filter out the first reads in this partition
        reads.zipPartitions(ctx.parallelize(firstReadInNextPartition, numParititions),
                (FlatMapFunction2<Iterator<GATKRead>, Iterator<List<GATKRead>>, GATKRead>) (it1, it2) -> {
                    PeekingIterator<GATKRead> current = Iterators.peekingIterator(it1);

                    final List<GATKRead> readsFromOtherPartition = new ArrayList<>();
                    it2.forEachRemaining(r -> {if (r!=null) readsFromOtherPartition.addAll(r);});
                    final int[] i = {0};

                    // test each read in this parition for equality with the other partition
                    current.forEachRemaining(currentRead -> {
                        for (GATKRead compareRead : readsFromOtherPartition) {
                            if (currentRead.equals(compareRead)) {
                                throw new GATKException(String.format("Read %s was duplicated in a partition at position &d in the partition", currentRead, i[0]));
                            }
                        }
                        i[0]++;
                    });
                    return null;
                });
    }
}