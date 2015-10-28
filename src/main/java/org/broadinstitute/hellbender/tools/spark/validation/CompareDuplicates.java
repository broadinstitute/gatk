package org.broadinstitute.hellbender.tools.spark.validation;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import htsjdk.samtools.SAMFileHeader;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.api.java.function.Function;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.SparkProgramGroup;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.engine.spark.datasources.ReadsSparkSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.transforms.markduplicates.MarkDuplicatesSparkUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadCoordinateComparator;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import scala.Tuple2;
import scala.Tuple3;

import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

@CommandLineProgramProperties(summary = "Compares two BAMs", oneLineSummary = "Compares two BAMs",
        programGroup = SparkProgramGroup.class)
public final class CompareDuplicates extends GATKSparkTool {
    private static final long serialVersionUID = 1L;

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



        JavaRDD<GATKRead> firstDupes = firstReads.filter(GATKRead::isDuplicate);
        JavaRDD<GATKRead> secondDupes = secondReads.filter(GATKRead::isDuplicate);
        long firstDupesCount = firstDupes.count();
        long secondDupesCount = secondDupes.count();
        System.out.println("first and second: " + firstDupesCount + "," + secondDupesCount);
        if (firstDupesCount != secondDupesCount) {
            throw new GATKException("BAMs have different number of total duplicates: " + firstDupesCount + "," + secondDupesCount);
        }
        JavaPairRDD<GATKRead, GATKRead> zip = firstReads.zip(secondReads);


        Broadcast<SAMFileHeader> broadcastHeader = ctx.broadcast(getHeaderForReads());

        JavaPairRDD<String, Iterable<GATKRead>> firstKeyed = firstReads.groupBy(new GroupKey());
        JavaPairRDD<String, Iterable<GATKRead>> secondKeyed = secondReads.groupBy(new GroupKey());

                JavaPairRDD < String, Tuple2 < Iterable < GATKRead >, Iterable < GATKRead >>> join = firstKeyed.join(secondKeyed);
        JavaPairRDD<String, Tuple3<String, Iterable<GATKRead>, Iterable<GATKRead>>> tagged = join.mapToPair(v1 -> {
            SAMFileHeader header = broadcastHeader.getValue();
            Iterable<GATKRead> firstReads1 = v1._2()._1();
            Iterable<GATKRead> secondReads1 = v1._2()._2();
            String str = Utils.getDupes(firstReads1, secondReads1, header); //new Dup().call(new Tuple3<>(firstReads1, secondReads1, header));
            return new Tuple2<>(v1._1(), new Tuple3<>(str, firstReads1, secondReads1));
        });

        Map<String, Integer> tagCountMap = tagged.mapToPair(v1 ->
                new Tuple2<>(v1._2()._1(), 1)).reduceByKey((v1, v2) -> v1 + v2).collectAsMap();
        Set<String> enums = Sets.newHashSet("EQUAL", "SIZE_UNEQUAL", "READ_MISMATCH", "DIFF_NUM_DUPES", "DUPES_MISMATCH");
        enums.forEach(s -> System.out.println(s + ": " + tagCountMap.getOrDefault(s, 0)));
        System.out.println("*************************");
        for (Map.Entry<String, Integer> e : tagCountMap.entrySet()) {
            System.out.println(e.getKey() + ": " + e.getValue());
        }

    }
}

class GroupKey implements Function<GATKRead, String> {
    private static final long serialVersionUID = 1L;
    @Override
    public String call(GATKRead v1) throws Exception {
        return v1.getContig() + "," + ReadUtils.getStrandedUnclippedStart(v1) + "," + ReadUtils.readHasMappedMate(v1);
    }
}


class Utils {
    static Tuple2<Integer, Integer> unpairedScoreMatches(
            List<Tuple2<Iterable<GATKRead>, Iterable<GATKRead>>> joinedUpairedReads) {
        // For unpaired, calculate the scores.
        int sameScore = 0;
        int differentScore = 0;
        for (Tuple2<Iterable<GATKRead>, Iterable<GATKRead>> entry : joinedUpairedReads) {
            // Find the score of the non-duplicate in the first BAM
            //System.out.println("********** Starting new location **********");
            int firstScore = -1;
            for (GATKRead first : entry._1()) {
                if (!first.isDuplicate()) {
                    if (firstScore == -1) {
                        firstScore = MarkDuplicatesSparkUtils.score(first);
                    } else {
                        System.out.println("...there should only be one non-duplicate, " + first);
                    }
                }
                //System.out.println("first: " + first + "," + first.isDuplicate() + "," + MarkDuplicatesSparkUtils.score(first));
            }
            int secondScore = -1;
            for (GATKRead second : entry._2()) {
                if (!second.isDuplicate()) {
                    if (secondScore == -1) {
                        secondScore = MarkDuplicatesSparkUtils.score(second);
                    } else {
                        System.out.println("...there should only be one non-duplicate, " + second);
                    }
                }
                //System.out.println("second: " + second + "," + second.isDuplicate() + "," + MarkDuplicatesSparkUtils.score(second));
            }
            if (firstScore == -1) {
                throw new GATKException("there should be at least one non-duplicate from the first bam");
            }
            if (secondScore == -1) {
                throw new GATKException("there should be at least one non-duplicate from the second bam, exiting...");
            }
            if (firstScore == secondScore) {
                sameScore++;
            } else {
                differentScore++;
            }
        }
        return new Tuple2<>(sameScore, differentScore);
    }

    static String getDupes(Iterable<GATKRead> f, Iterable<GATKRead> s, SAMFileHeader header) {
        List<GATKRead> first = Lists.newArrayList(f);
        List<GATKRead> second = Lists.newArrayList(s);
        if (first.size() != second.size()) {
            return  "SIZE_UNEQUAL";
        }
        int size = first.size();
        if (size > 1) {
            first.sort(new ReadCoordinateComparator(header));
            second.sort(new ReadCoordinateComparator(header));
        }

        HashSet<GATKRead> firstDupes = Sets.newHashSet();
        HashSet<GATKRead> secondDupes = Sets.newHashSet();
        for (int i = 0; i < size; ++i) {
            GATKRead firstRead = first.get(i);
            GATKRead secondRead = second.get(i);
            if (!firstRead.getName().equals(secondRead.getName())) {
                return "READ_MISMATCH";
            }
            if (firstRead.isDuplicate()) {
                firstDupes.add(firstRead);
            }
            if (secondRead.isDuplicate()) {
                secondDupes.add(secondRead);
            }
        }
        if (firstDupes.size() != secondDupes.size()) {
            return "DIFF_NUM_DUPES";
        }
        if (firstDupes.isEmpty()) {
            return "EQUAL";
        }
        if (!firstDupes.containsAll(secondDupes)) {
            return "DUPES_MISMATCH";
        }
        return "EQUAL";
    }
}

/*
class Dup implements Function<Tuple3<Iterable<GATKRead>, Iterable<GATKRead>, SAMFileHeader>, String> {
    @Override
    public String call(Tuple3<Iterable<GATKRead>, Iterable<GATKRead>, SAMFileHeader> v1) throws Exception {
        SAMFileHeader header = v1._3();
        List<GATKRead> first = Lists.newArrayList(v1._1());
        List<GATKRead> second = Lists.newArrayList(v1._2());
        if (first.size() != second.size()) {
            return  "SIZE_UNEQUAL";
        }
        int size = first.size();
        if (size > 1) {
            first.sort(new ReadCoordinateComparator(header));
            second.sort(new ReadCoordinateComparator(header));
        }

        HashSet<GATKRead> firstDupes = Sets.newHashSet();
        HashSet<GATKRead> secondDupes = Sets.newHashSet();
        for (int i = 0; i < size; ++i) {
            GATKRead firstRead = first.get(i);
            GATKRead secondRead = second.get(i);
            if (!firstRead.getName().equals(secondRead.getName())) {
                return "READ_MISMATCH";
            }
            if (firstRead.isDuplicate()) {
                firstDupes.add(firstRead);
            }
            if (secondRead.isDuplicate()) {
                secondDupes.add(secondRead);
            }
        }
        if (firstDupes.size() != secondDupes.size()) {
            return "DIFF_NUM_DUPES";
        }
        if (firstDupes.isEmpty()) {
            return "EQUAL";
        }
        if (!firstDupes.containsAll(secondDupes)) {
            return "DUPES_MISMATCH";
        }
        return "EQUAL";
    }
}*/

