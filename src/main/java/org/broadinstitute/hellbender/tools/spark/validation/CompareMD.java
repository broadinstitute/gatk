package org.broadinstitute.hellbender.tools.spark.validation;

import com.google.common.collect.ListMultimap;
import com.google.common.collect.Lists;
import com.google.common.collect.Multimaps;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.api.java.function.Function;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.SparkProgramGroup;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.engine.spark.datasources.ReadsSparkSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.dataflow.transforms.markduplicates.PairedEnds;
import org.broadinstitute.hellbender.tools.dataflow.transforms.markduplicates.ReadsKey;
import org.broadinstitute.hellbender.tools.spark.transforms.markduplicates.MarkDuplicatesSparkUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import scala.Tuple2;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

@CommandLineProgramProperties(summary = "Compares two BAMs", oneLineSummary = "Compares two BAMs",
        programGroup = SparkProgramGroup.class)

public final class CompareMD extends GATKSparkTool {
    private static final long serialVersionUID = 1L;

    private static final String FirstDuplicate = "FD";
    private static final String SecondDuplicate = "SD";

    @Override
    public boolean requiresReads() { return true; }

    @Argument(doc="the second BAM", shortName = "I2", fullName = "input2", optional = false)
    protected String input2;

    @Override
    protected void runTool(final JavaSparkContext ctx) {
        JavaRDD<GATKRead> firstReads = getReads().map(v1 -> {
            v1.clearAttributes();
            return v1;
        });
        ReadsSparkSource readsSource2 = new ReadsSparkSource(ctx);
        JavaRDD<GATKRead> secondReads = readsSource2.getParallelReads(input2, getIntervals()).map(v1 -> {
            v1.clearAttributes();
            return v1;
        });


        long firstBamSize = firstReads.count();
        long secondBamSize = secondReads.count();
        if (firstBamSize != secondBamSize) {
            throw new UserException("input bams have different numbers of total reads: "
                    + firstBamSize + "," + secondBamSize);
        }
        System.out.println("processing bams with " + firstBamSize + " reads");

        JavaRDD<GATKRead> firstPrimary = firstReads.filter(v1 -> !isNonPrimary(v1));
        JavaRDD<GATKRead> secondPrimary = secondReads.filter(v1 -> !isNonPrimary(v1));
        // TODO: Verify all non-primary are non-duplicates.
        long firstNonPrimaryDupes = firstReads.filter(v1 -> isNonPrimary(v1) && v1.isDuplicate()).count();
        long secondNonPrimaryDupes = secondReads.filter(v1 -> isNonPrimary(v1) && v1.isDuplicate()).count();
        System.out.println("firstNonPrimaryDupes, secondNonPrimaryDupes: " + firstNonPrimaryDupes + "," + secondNonPrimaryDupes);

        // Group the reads together and add tags for first and second bam duplicates
        JavaPairRDD < String, GATKRead > firstNamed = firstPrimary.mapToPair(read -> new Tuple2<>(read.toString(), read));
        firstNamed.groupByKey().foreach(stringIterableTuple2 -> {
            assert Lists.newArrayList(stringIterableTuple2._2()).size() == 1;
        });
        JavaPairRDD<String, GATKRead> secondNamed = secondPrimary.mapToPair(read -> new Tuple2<>(read.toString(), read));
        secondNamed.groupByKey().foreach(stringIterableTuple2 -> {
            assert Lists.newArrayList(stringIterableTuple2._2()).size() == 1;
        });
        JavaPairRDD<String, Tuple2<GATKRead, GATKRead>> joined = firstNamed.join(secondNamed);
        JavaRDD<GATKRead> readWithAttributes = joined.map(v1 -> {
            GATKRead firstRead = v1._2()._1();
            GATKRead secondRead = v1._2()._2();
            if (firstRead.isDuplicate()) {
                firstRead.setAttribute(FirstDuplicate, "");
            }
            if (secondRead.isDuplicate()) {
                firstRead.setAttribute(SecondDuplicate, "");
            }
            return firstRead;
        });

        long fd = readWithAttributes.filter(v1 -> v1.hasAttribute(FirstDuplicate)).count();
        long mismatches = readWithAttributes.filter(v1 -> {
            if (v1.hasAttribute(FirstDuplicate) && v1.hasAttribute(SecondDuplicate)) {
                return false;
            }
            if (!v1.hasAttribute(FirstDuplicate) && !v1.hasAttribute(SecondDuplicate)) {
                return false;
            }
            return true;
        }).count();
        System.out.println("fd, mismatches: " + fd + ", " + mismatches);
        int parallelism = readWithAttributes.partitions().size();

        // Reproduce the grouping of MarkDuplicates to find how the duplicates differ
        JavaPairRDD<String, Iterable<PairedEnds>> groupedPairs = MarkDuplicatesSparkUtils.groupReads(getHeaderForReads(), readWithAttributes, parallelism);

        JavaRDD<String> stats = groupedPairs.map(v1 -> {
            Iterable<PairedEnds> pairedEnds = v1._2();


            // Each key corresponds to either fragments or paired ends, not a mixture of both.
            if (ReadsKey.isFragment(v1._1())) { // fragments
                return compareDuplicateFragments(pairedEnds);
            }
            return comparePairedReads(pairedEnds);
        });
        Map<String, Integer> Tags = stats.mapToPair(s -> new Tuple2<>(s, 1)).reduceByKey((v1, v2) -> v1 + v2).collectAsMap();
        for (Map.Entry<String, Integer> entry : Tags.entrySet()) {
            System.out.println(entry.getKey() + "," + entry.getValue());
        }

    }
    private static boolean isNonPrimary(GATKRead read) {
        return read.isSecondaryAlignment() || read.isSupplementaryAlignment() || read.isUnmapped();
    }

    public static String compareDuplicateFragments(Iterable<PairedEnds> pairedEnds) {
        List<PairedEnds> pairs = Lists.newArrayList(pairedEnds);
        List<PairedEnds> fragments = pairs.stream().filter(pe -> !ReadUtils.readHasMappedMate(pe.first())).collect(Collectors.toList());
        if (!sameNumberDupes(fragments)) {
            return "DIFFERENT_NUM_DUPES__FRAGMENTS";
        }
        return getLabel(fragments) + "__FRAGMENTS";
    }

    public static Boolean sameNumberDupes(final List<PairedEnds> pairs) {
        int firstBamDupes = 0;
        int secondBamDupes = 0;

        for (final PairedEnds pair : pairs) {
            GATKRead r = pair.first();
            if (r.hasAttribute(FirstDuplicate)) {
                firstBamDupes++;
            }
            if (r.hasAttribute(SecondDuplicate)) {
                secondBamDupes++;
            }
            if (pair.second() != null) {
                GATKRead r2 = pair.second();
                if (r2.hasAttribute(FirstDuplicate)) {
                    firstBamDupes++;
                }
                if (r2.hasAttribute(SecondDuplicate)) {
                    secondBamDupes++;
                }
            }
        }
        return firstBamDupes == secondBamDupes;
    }

    public static String getLabel(final List<PairedEnds> pairs) {
        // See if all pairs are marked identically.
        Boolean identical = true;
        for (final PairedEnds p : pairs) {
            //for (GATKRead first : firstReads) {
            GATKRead first = p.first();
            if (first.hasAttribute(FirstDuplicate) && !first.hasAttribute(SecondDuplicate)) {
                identical = false;
            }
            if (!first.hasAttribute(FirstDuplicate) && first.hasAttribute(SecondDuplicate)) {
                identical = false;
            }
            if (p.second() != null) {
                GATKRead second = p.second();
                if (second.hasAttribute(FirstDuplicate) && !second.hasAttribute(SecondDuplicate)) {
                    identical = false;
                }
                if (!second.hasAttribute(FirstDuplicate) && second.hasAttribute(SecondDuplicate)) {
                    identical = false;
                }
            }
        }
        if (identical) {
            return "NO_DIFF";
        }

        // There are the number of duplicates, but not everything is marked identically.
        // Check to see if the picked reads have the same score (implying the difference is from tie-breaking).

        int firstReadScore = -1;
        int secondReadScore = -1;
        for (final PairedEnds p : pairs) {
            if (!p.first().hasAttribute(FirstDuplicate)) {
                if (firstReadScore != -1) {
                    throw new GATKException("two non-duplicates reads found (first BAM)");
                }
                if (p.second() != null) {
                    firstReadScore = p.score();
                } else {
                    firstReadScore = MarkDuplicatesSparkUtils.score(p.first());
                }
            }
            if (!p.first().hasAttribute(SecondDuplicate)) {
                if (secondReadScore != -1) {
                    throw new GATKException("two non-duplicates reads found (second BAM)");
                }
                if (p.second() != null) {
                    secondReadScore = p.score();
                } else {
                    secondReadScore = MarkDuplicatesSparkUtils.score(p.first());
                }
            }
        }
        if (firstReadScore != secondReadScore) {
            return "DIFFERENT_TOP_SCORES";
        }
        return "SAME_TOP_SCORES";
    }

    public static String comparePairedReads(Iterable<PairedEnds> pairedEnds) {

        final ListMultimap<Boolean, PairedEnds> paired = Multimaps.index(pairedEnds, pair -> pair.second() != null);
        List<PairedEnds> singles = paired.get(false);
        List<PairedEnds> pairs = paired.get(true);

        // First, count the total duplicates for this key
        if (!sameNumberDupes(singles) || !sameNumberDupes(pairs)) {
            return "DIFFERENT_NUM_DUPES__PAIRED";
        }


        for (final PairedEnds pair : singles) {
            if (pair.first().hasAttribute(FirstDuplicate) || pair.first().hasAttribute(SecondDuplicate)) {
                // We expect unpaired reads (note: not fragments), to be left alone
                return "UNPAIRED_DUPLICATE_MARKED__PAIRED";
            }
        }

        return getLabel(pairs) + "__PAIRED";
    }
}
