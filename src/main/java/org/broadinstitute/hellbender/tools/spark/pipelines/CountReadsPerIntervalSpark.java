package org.broadinstitute.hellbender.tools.spark.pipelines;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.OverlapDetector;
import org.apache.spark.api.java.JavaPairRDD;
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
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.spark.SparkUtils;
import scala.Tuple2;

import java.io.PrintStream;
import java.util.*;
import java.util.stream.Collectors;

@CommandLineProgramProperties(summary = "Counts reads per interval in the input SAM/BAM",
        oneLineSummary = "CountReads per interval on Spark",
        programGroup = SparkProgramGroup.class)
public final class CountReadsPerIntervalSpark extends GATKSparkTool {

    private static final long serialVersionUID = 1L;

    @Override
    public boolean requiresReads() { return true; }

    @Argument(doc = "uri for the output file: a local file path",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            optional = false)
    public String out;

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
        chromosomeSizes.put("1", 248956422);
        chromosomeSizes.put("2", 242193529);
        chromosomeSizes.put("3", 198295559);
        chromosomeSizes.put("4", 190214555);
        chromosomeSizes.put("5", 181538259);
        chromosomeSizes.put("6", 170805979);
        chromosomeSizes.put("7", 159345973);
        chromosomeSizes.put("8", 145138636);
        chromosomeSizes.put("9", 138394717);
        chromosomeSizes.put("10", 133797422);
        chromosomeSizes.put("11", 135086622);
        chromosomeSizes.put("12", 133275309);
        chromosomeSizes.put("13", 114364328);
        chromosomeSizes.put("14", 107043718);
        chromosomeSizes.put("15", 101991189);
        chromosomeSizes.put("16", 90338345);
        chromosomeSizes.put("17", 83257441);
        chromosomeSizes.put("18", 80373285);
        chromosomeSizes.put("19", 58617616);
        chromosomeSizes.put("20", 64444167);
        chromosomeSizes.put("21", 46709983);
        chromosomeSizes.put("22", 50818468);
        chromosomeSizes.put("X", 156040895);
        chromosomeSizes.put("Y", 57227415);
        chromosomeSizes.put("M", 16569);
        for (Map.Entry<String, Integer> e : chromosomeSizes.entrySet()) {
            String contig = e.getKey();
            int size = e.getValue();
            for (int i = 0; i <= size/(regionSize * regionRepeat); i++) {
                int start = i * regionSize * regionRepeat + 1;
                int end = start + regionSize;
                intervals.add(new SimpleInterval(contig, start, end));
            }
        }

        JavaPairRDD<Locatable, Integer> intervalCountsRdd;

        if (shuffle) {
            intervalCountsRdd = SparkUtils.joinOverlappingShuffle(ctx, mappedReads, GATKRead.class, intervals,
                    new CountOverlappingReadsFlatFunction()).mapToPair(t -> t);
        } else {
            intervalCountsRdd = SparkUtils.joinOverlapping(ctx, mappedReads, GATKRead.class, intervals,
                    new CountOverlappingReadsFunction()).mapToPair(t -> t);
        }

        List<Tuple2<Locatable, Integer>> collect = new ArrayList<>(intervalCountsRdd.collect());
        Collections.sort(collect, (o1, o2) -> {
            Locatable l1 = o1._1();
            Locatable l2 = o2._1();
            int c = l1.getContig().compareTo(l2.getContig());
            if (c != 0) {
                return c;
            }
            c = Integer.compare(l1.getStart(), l2.getStart());
            if (c != 0) {
                return c;
            }
            return Integer.compare(l1.getEnd(), l2.getEnd());
        });
        //collect.forEach(t -> System.out.println(t));

        if(out != null) {
            try (final PrintStream ps = new PrintStream(BucketUtils.createFile(out, getAuthenticatedGCSOptions()))) {
                collect.forEach(t -> ps.println(t));
            }
        }
    }

    @SuppressWarnings("unchecked")
    private static class CountOverlappingReadsFunction implements FlatMapFunction2<Iterator<GATKRead>, Iterator<Locatable>, Tuple2<Locatable, Integer>> {
        private static final long serialVersionUID = 1L;

        @Override
        public Iterable<Tuple2<Locatable, Integer>> call(Iterator<GATKRead> readIterator,
                                                         Iterator<Locatable> intervalIterator) {
            List<Locatable> intervals = Lists.newArrayList(intervalIterator);
            OverlapDetector<Locatable> overlapDetector = OverlapDetector.create(intervals);

            Map<Locatable, Integer> counts = Maps.newLinkedHashMap();
            while (readIterator.hasNext()) {
                GATKRead read = readIterator.next();
                Set<Locatable> overlaps = overlapDetector.getOverlaps(read);
                for (Locatable overlap : overlaps) {
                    Integer count = counts.get(overlap);
                    counts.put(overlap, count == null ? 1 : count + 1);
                }
            }
            return counts.entrySet().stream()
                    .map(entry -> new Tuple2<>(entry.getKey(), entry.getValue()))
                    .collect(Collectors.toList());
        }
    }

    private static class CountOverlappingReadsFlatFunction implements FlatMapFunction<Tuple2<Locatable, Iterable<GATKRead>>, Tuple2<Locatable, Integer>> {
        private static final long serialVersionUID = 1L;

        @Override
        public Iterable<Tuple2<Locatable, Integer>> call(Tuple2<Locatable, Iterable<GATKRead>> t) throws Exception {
            Locatable interval = t._1;
            int count = 0;
            OverlapDetector<Locatable> overlapDetector = OverlapDetector.create(ImmutableList.of(interval));
            for (GATKRead read : t._2) {
                count += overlapDetector.getOverlaps(read).size();
            }
            return ImmutableList.of(new Tuple2<>(interval, count));
        }
    }

}
