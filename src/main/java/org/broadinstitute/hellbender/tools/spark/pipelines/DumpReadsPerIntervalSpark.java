package org.broadinstitute.hellbender.tools.spark.pipelines;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Iterables;
import com.google.common.collect.Iterators;
import com.google.common.collect.Lists;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.OverlapDetector;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.api.java.function.FlatMapFunction;
import org.apache.spark.api.java.function.FlatMapFunction2;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.SparkProgramGroup;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.spark.SparkUtils;
import scala.Tuple2;

import java.util.*;

@CommandLineProgramProperties(summary = "Counts reads per interval in the input SAM/BAM",
        oneLineSummary = "CountReads per interval on Spark",
        programGroup = SparkProgramGroup.class)
public final class DumpReadsPerIntervalSpark extends GATKSparkTool {

    private static final long serialVersionUID = 1L;

    @Override
    public boolean requiresReads() { return true; }

    @Argument(doc="the output file path", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, optional = false)
    protected String outputFile;

    @Argument(doc = "whether to use the shuffle implementation or not", shortName = "shuffle", fullName = "shuffle", optional = true)
    public boolean shuffle = false;

    @Override
    protected void runTool(final JavaSparkContext ctx) {
        final JavaRDD<GATKRead> reads = getReads();
        JavaRDD<GATKRead> mappedReads = reads.filter(read -> !read.isUnmapped());

        final List<Locatable> intervals = new ArrayList<>();
        int regionSize = 4000;
        int regionRepeat = 3; // fraction of regions 3 = 1/3 etc
        for (SAMSequenceRecord sequenceRecord : getBestAvailableSequenceDictionary().getSequences()) {
            String contig = sequenceRecord.getSequenceName();
            int size = sequenceRecord.getSequenceLength();
            for (int i = 0; i <= size/(regionSize * regionRepeat); i++) {
                int start = i * regionSize * regionRepeat + 1;
                int end = start + regionSize;
                intervals.add(new SimpleInterval(contig, start, end));
            }
        }

        JavaRDD<GATKRead> collectedReads;

        if (shuffle) {
            collectedReads = SparkUtils.joinOverlappingShuffle(ctx, mappedReads, GATKRead.class, intervals,
                    new DumpOverlappingReadsFlatFunction());
        } else {
            collectedReads = SparkUtils.joinOverlapping(ctx, mappedReads, GATKRead.class, getBestAvailableSequenceDictionary(),
                    intervals, new DumpOverlappingReadsFunction());
        }
        shardedOutput = true;
        writeReads(ctx, outputFile, collectedReads);
    }

    @SuppressWarnings("unchecked")
    private static class DumpOverlappingReadsFunction implements FlatMapFunction2<Iterator<GATKRead>, Iterator<Locatable>, GATKRead> {
        private static final long serialVersionUID = 1L;

        @Override
        public Iterable<GATKRead> call(Iterator<GATKRead> readIterator,
                                                         Iterator<Locatable> intervalIterator) {
            List<Locatable> intervals = Lists.newArrayList(intervalIterator);
            OverlapDetector<Locatable> overlapDetector = OverlapDetector.create(intervals);
            return () -> Iterators.filter(readIterator, overlapDetector::overlapsAny);
        }
    }

    private static class DumpOverlappingReadsFlatFunction implements FlatMapFunction<Tuple2<Locatable, Iterable<GATKRead>>, GATKRead> {
        private static final long serialVersionUID = 1L;

        @Override
        public Iterable<GATKRead> call(Tuple2<Locatable, Iterable<GATKRead>> t) throws Exception {
            Locatable interval = t._1;
            OverlapDetector<Locatable> overlapDetector = OverlapDetector.create(ImmutableList.of(interval));
            return Iterables.filter(t._2(), overlapDetector::overlapsAny);
        }
    }

}
