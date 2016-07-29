package org.broadinstitute.hellbender.tools.spark.pipelines;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Iterables;
import com.google.common.collect.Iterators;
import com.google.common.collect.Lists;
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

import javax.annotation.Nullable;
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
        Map<String, Integer> chromosomeSizes = new LinkedHashMap<>();
        int regionSize = 4000;
        int regionRepeat = 3; // fraction of regions 3 = 1/3 etc
        chromosomeSizes.put("chr1", 248956422);
        chromosomeSizes.put("chr2", 242193529);
        chromosomeSizes.put("chr3", 198295559);
        chromosomeSizes.put("chr4", 190214555);
        chromosomeSizes.put("chr5", 181538259);
        chromosomeSizes.put("chr6", 170805979);
        chromosomeSizes.put("chr7", 159345973);
        chromosomeSizes.put("chr8", 145138636);
        chromosomeSizes.put("chr9", 138394717);
        chromosomeSizes.put("chr10", 133797422);
        chromosomeSizes.put("chr11", 135086622);
        chromosomeSizes.put("chr12", 133275309);
        chromosomeSizes.put("chr13", 114364328);
        chromosomeSizes.put("chr14", 107043718);
        chromosomeSizes.put("chr15", 101991189);
        chromosomeSizes.put("chr16", 90338345);
        chromosomeSizes.put("chr17", 83257441);
        chromosomeSizes.put("chr18", 80373285);
        chromosomeSizes.put("chr19", 58617616);
        chromosomeSizes.put("chr20", 64444167);
        chromosomeSizes.put("chr21", 46709983);
        chromosomeSizes.put("chr22", 50818468);
        chromosomeSizes.put("chrX", 156040895);
        chromosomeSizes.put("chrY", 57227415);
        chromosomeSizes.put("chrM", 16569);
        for (Map.Entry<String, Integer> e : chromosomeSizes.entrySet()) {
            String contig = e.getKey();
            int size = e.getValue();
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
            collectedReads = SparkUtils.joinOverlapping(ctx, mappedReads, GATKRead.class, intervals,
                    new DumpOverlappingReadsFunction());
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
