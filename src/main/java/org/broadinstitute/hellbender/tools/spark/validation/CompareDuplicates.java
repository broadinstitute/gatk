package org.broadinstitute.hellbender.tools.spark.validation;

import com.google.common.collect.Iterables;
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
        JavaPairRDD<String, Tuple3<String, Iterable<GATKRead>, Iterable<GATKRead>>> equal =
                tagged.filter(v1 -> v1._2()._1().equals("EQUAL"));
        JavaPairRDD<String, Tuple3<String, Iterable<GATKRead>, Iterable<GATKRead>>> unequal =
                tagged.filter(v1 -> !v1._2()._1().equals("EQUAL"));

        long diffNumDupes = unequal.filter(v1 -> v1._2()._1().equals("DIFF_NUM_DUPES")).count();

        Map<String, Tuple3<String, Iterable<GATKRead>, Iterable<GATKRead>>> dupesMismatch = unequal.filter(v1 -> v1._2()._1().equals("DUPES_MISMATCH")).collectAsMap();
        long dupesMismatchSize = dupesMismatch.size();

        // For all of the mismatch cases, split them into paired and unpaired
        List<Tuple2<Iterable<GATKRead>, Iterable<GATKRead>>> joinedUpairedReads = Lists.newArrayList();
        List<Tuple2<Iterable<GATKRead>, Iterable<GATKRead>>> joinedPairedReads = Lists.newArrayList();
        for (Map.Entry<String, Tuple3<String, Iterable<GATKRead>, Iterable<GATKRead>>> entry : dupesMismatch.entrySet()) {
            Tuple3<String, Iterable<GATKRead>, Iterable<GATKRead>> item = entry.getValue();
            GATKRead first = Iterables.getFirst(item._2(), null); // the item must be non-null
            if (ReadUtils.readHasMappedMate(first)) {
                joinedPairedReads.add(new Tuple2<>(item._2(), item._3()));
            } else {
                joinedUpairedReads.add(new Tuple2<>(item._2(), item._3()));

            }
        }

        // For the unpaired reads, see if the scores match.
        Tuple2<Integer, Integer> unpairedSameDifferentScore = Utils.unpairedScoreMatches(joinedUpairedReads);
        int sameScore = unpairedSameDifferentScore._1();
        int differentScore = unpairedSameDifferentScore._2();


        // For paired reads, we need to find the mate (so it too can be score).
        Set<String> firstDuplicates = Sets.newHashSet();
        Set<String> secondDuplicates = Sets.newHashSet();
        for (Tuple2<Iterable<GATKRead>, Iterable<GATKRead>> entry : joinedPairedReads) {
            Iterable<GATKRead> iFirst = entry._1();
            Iterable<GATKRead> iSecond = entry._2();
            for (GATKRead r : iFirst) {
                firstDuplicates.add(r.getName());
            }
            for (GATKRead r : iSecond) {
                secondDuplicates.add(r.getName());
            }
        }
        Broadcast<Set<String>> firstDupesBroadcast = ctx.broadcast(firstDuplicates);
        Broadcast<Set<String>> secondDupesBroadcast = ctx.broadcast(secondDuplicates);
        Map<String, Iterable<GATKRead>> firstDupPairsMap = firstReads.filter(v1 -> {
            Set<String> firstDupes = firstDupesBroadcast.getValue();
            return firstDupes.contains(v1.getName());
        }).groupBy(GATKRead::getName).collectAsMap();



        Map<String, Iterable<GATKRead>> secondDupPairsMap = secondReads.filter(v1 -> {
            Set<String> secondDupes = secondDupesBroadcast.getValue();
            return secondDupes.contains(v1.getName());
        }).groupBy(GATKRead::getName).collectAsMap();


        int samePairedScore = 0;
        int differentPairedScore = 0;
        // Look at paired dupes scores
        for (Tuple2<Iterable<GATKRead>, Iterable<GATKRead>> entry : joinedPairedReads) {
            int firstScore = -1;
            for (GATKRead first : entry._1()) {
                if (!first.isDuplicate()) {
                    if (firstScore == -1) {
                        Iterable<GATKRead> pair = firstDupPairsMap.get(first.getName());
                        List<GATKRead> lPair = Lists.newArrayList(pair);
                        if (lPair.size() != 2) {
                            throw new GATKException("expected to find two reads from mate map, found " + lPair.size());
                        }
                        firstScore = MarkDuplicatesSparkUtils.score(lPair.get(0)) + MarkDuplicatesSparkUtils.score(lPair.get(1));
                    } else {
                        System.out.println("...there should only be one paired non-duplicate (1), " + first);
                        for (GATKRead f : entry._1()) {
                            System.out.println("read: " + f + ", " + f.isDuplicate());
                        }
                    }
                }
            }
            int secondScore = -1;
            for (GATKRead second : entry._2()) {
                if (!second.isDuplicate()) {
                    if (secondScore == -1) {
                        Iterable<GATKRead> pair = secondDupPairsMap.get(second.getName());
                        List<GATKRead> lPair = Lists.newArrayList(pair);
                        if (lPair.size() != 2) {
                            throw new GATKException("expected to find two reads from mate map, found " + lPair.size());
                        }
                        secondScore = MarkDuplicatesSparkUtils.score(lPair.get(0)) + MarkDuplicatesSparkUtils.score(lPair.get(1));
                    } else {
                        System.out.println("...there should only be one paired non-duplicate (2), " + second);
                    }
                }
            }

            if (firstScore == -1) {
                throw new GATKException("there should be at least one non-duplicate from the first bam");
            }
            if (secondScore == -1) {
                throw new GATKException("there should be at least one non-duplicate from the second bam, exiting...");
            }
            if (firstScore == secondScore) {
                samePairedScore++;
            } else {
                differentPairedScore++;
            }
        }



        JavaRDD < GATKRead > firstDupes = firstReads.filter(GATKRead::isDuplicate);
        JavaRDD<GATKRead> secondDupes = secondReads.filter(GATKRead::isDuplicate);

        /*
        JavaRDD<String> firstDistinct = firstDupes.map((Function<GATKRead, String>) (read) -> {
            boolean paired = read.isPaired();
            boolean mateIsReverseStrand = false;
            if (paired) {
                mateIsReverseStrand = read.mateIsReverseStrand();
            }
            return read.getName() + "," + paired + "," + mateIsReverseStrand;
        }).distinct();
        JavaRDD<String> secondDistinct = secondDupes.map((Function<GATKRead, String>) (read) -> {
            boolean paired = read.isPaired();
            boolean mateIsReverseStrand = false;
            if (paired) {
                mateIsReverseStrand = read.mateIsReverseStrand();
            }
            return read.getName() + "," + paired + "," + mateIsReverseStrand;
        }).distinct();*/
        long firstDupesCount = firstDupes.count();
        long secondDupesCount = secondDupes.count();
        System.out.println("first and second: " + firstDupesCount + "," + secondDupesCount);
        //System.out.println("first and second, total distinct: " + firstDistinct.count() + "," +
        //        secondDistinct.count() + "," + allDupes);
        System.out.println("no diff,diffs,diffNumDupes,dupesMismatch: " + equal.count() + "," + unequal.count() + "," +
                diffNumDupes + "," + dupesMismatchSize);


        System.out.println("same score, different score: " + sameScore + "," + differentScore);
        System.out.println("same paired score, different paired score: " + samePairedScore + "," + differentPairedScore);

        /*
        for (Map.Entry<String, Tuple3<String, Iterable<GATKRead>, Iterable<GATKRead>>> entry : diffNumDupesMap.entrySet()) {
            System.out.println("**************************************");
            System.out.println("*************** " + entry.getKey() + "*************************");
            System.out.println("**************************************");
            List<GATKRead> iFirst = Lists.newArrayList(entry.getValue()._2());
            List<GATKRead> iSecond = Lists.newArrayList(entry.getValue()._3());
            System.out.println("first");
            for (GATKRead r : iFirst) {
                System.out.println(r.isDuplicate() + "," + ((SAMRecordToGATKReadAdapter) r).getSamRecord().getSAMString());
            }
            System.out.println("second");
            for (GATKRead r : iSecond) {
                System.out.println(r.isDuplicate() + "," + ((SAMRecordToGATKReadAdapter) r).getSamRecord().getSAMString());
            }
        }*/
        

        /*
        unequal.foreach(v1 -> System.out.println(v1._2()._1() + "," +
                Lists.newArrayList(v1._2()._2()).size() + "," +
                Lists.newArrayList(v1._2()._3()).size()));
        */
    }
}

class GroupKey implements Function<GATKRead, String> {

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

