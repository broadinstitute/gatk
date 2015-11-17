package org.broadinstitute.hellbender.tools.spark.validation;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import htsjdk.samtools.SAMFileHeader;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.api.java.function.Function;
import org.apache.spark.api.java.function.Function2;
import org.apache.spark.api.java.function.PairFunction;
import org.apache.spark.broadcast.Broadcast;
import org.apache.spark.storage.StorageLevel;
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

import java.io.Serializable;
import java.util.*;
import java.util.function.Consumer;

@CommandLineProgramProperties(summary = "Compares two BAMs", oneLineSummary = "Compares two BAMs",
        programGroup = SparkProgramGroup.class)
public final class CompareDuplicates extends GATKSparkTool {
    private static final long serialVersionUID = 1L;

    private final double sampleRate = 0.01;
    private final int pathologicalThreshold = 2;//100000;
    @Override
    public boolean requiresReads() { return true; }

    @Argument(doc="the second BAM", shortName = "I2", fullName = "input2", optional = false)
    protected String input2;


    @Override
    protected void runTool(final JavaSparkContext ctx) {
        JavaRDD<GATKRead> firstReads = getReads().filter(v1 -> {
            if (Utils.isNonPrimary(v1) && v1.isDuplicate()) {
                throw new RuntimeException("found a non-primary read marked as a duplicate in the first bam (this shouldn't happen)");
            }
            return !Utils.isNonPrimary(v1);
        });
        firstReads.foreach(GATKRead::clearAttributes);
        //firstReads.persist(StorageLevel.MEMORY_AND_DISK());
        ReadsSparkSource readsSource2 = new ReadsSparkSource(ctx);
        JavaRDD<GATKRead> secondReads = readsSource2.getParallelReads(input2, getIntervals()).filter(v1 -> {
            if (Utils.isNonPrimary(v1) && v1.isDuplicate()) {
                throw new RuntimeException("found a non-primary read marked as a duplicate in the second bam (this shouldn't happen)");
            }
            return !Utils.isNonPrimary(v1);
        });
        secondReads.foreach(GATKRead::clearAttributes);
        //secondReads.persist(StorageLevel.MEMORY_AND_DISK());

        
        long firstBamSize = firstReads.count();
        JavaRDD<GATKRead> firstDupes = firstReads.filter(GATKRead::isDuplicate);
        long firstDupesCount = firstDupes.count();

        long secondBamSize = secondReads.count();
        JavaRDD<GATKRead> secondDupes = secondReads.filter(GATKRead::isDuplicate);
        long secondDupesCount = secondDupes.count();

        if (firstBamSize != secondBamSize) {
            throw new UserException("input bams have different numbers of total reads: "
                    + firstBamSize + "," + secondBamSize);
        }
        System.out.println("processing bams with " + firstBamSize + " mapped reads");


        System.out.println("first and second: " + firstDupesCount + "," + secondDupesCount);
        if (firstDupesCount != secondDupesCount) {
            throw new GATKException("BAMs have different number of total duplicates: " + firstDupesCount + "," + secondDupesCount);
        }

        JavaPairRDD<String, GATKRead> firstKeyed = firstReads.mapToPair(read -> new Tuple2<>(Utils.groupKey(read), read));
        JavaPairRDD<String, GATKRead> secondKeyed = secondReads.mapToPair(read -> new Tuple2<>(Utils.groupKey(read), read));
        JavaPairRDD<String, Tuple2<Iterable<GATKRead>, Iterable<GATKRead>>> cogroup = firstKeyed.cogroup(secondKeyed);

        Broadcast<SAMFileHeader> broadcastHeader = ctx.broadcast(getHeaderForReads());
        JavaRDD<String> tagged = cogroup.map(v1 -> {
            SAMFileHeader header = broadcastHeader.getValue();

            Iterable<GATKRead> firstReads1 = v1._2()._1();
            Iterable<GATKRead> secondReads1 = v1._2()._2();

            return Utils.getDupes(firstReads1, secondReads1, header);
        });

        Map<String, Integer> tagCountMap = tagged.mapToPair(v1 ->
                new Tuple2<>(v1, 1)).reduceByKey((v1, v2) -> v1 + v2).collectAsMap();
        Set<String> enums = Sets.newHashSet("EQUAL", "SIZE_UNEQUAL", "READ_MISMATCH", "DIFF_NUM_DUPES", "DUPES_MISMATCH");
        enums.forEach(s -> System.out.println(s + ": " + tagCountMap.getOrDefault(s, 0)));
        System.out.println("*************************");
        for (Map.Entry<String, Integer> e : tagCountMap.entrySet()) {
            System.out.println(e.getKey() + ": " + e.getValue());
        }

        try {
            Thread.sleep(120*1000);                 //1000 milliseconds is one second.
        } catch(InterruptedException ex) {
            Thread.currentThread().interrupt();
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


class Utils implements Serializable {
    private static final long serialVersionUID = 1L;

    static boolean isNonPrimary(GATKRead read) {
        return read.isSecondaryAlignment() || read.isSupplementaryAlignment() || read.isUnmapped();
    }

    static String groupKey(GATKRead r) {
        return r.getContig() + "," + ReadUtils.getStrandedUnclippedStart(r) + "," + ReadUtils.readHasMappedMate(r);
    }

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

    static String getDupes2(Iterable<GATKRead> f) {
        List<GATKRead> first = Lists.newArrayList(f);

        HashSet<GATKRead> firstDupes = Sets.newHashSet();
        HashSet<GATKRead> secondDupes = Sets.newHashSet();
        for (GATKRead r : first) {
            if (r.hasAttribute("FD")) {
                firstDupes.add(r);
            }
            if (r.hasAttribute("SD")) {
                secondDupes.add(r);
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

    void getStats(JavaRDD<GATKRead> toCache, double sampleRate) {
        toCache.persist(StorageLevel.MEMORY_AND_DISK());

        JavaRDD<GATKRead> sample = toCache.sample(false, sampleRate);
        JavaPairRDD<String, Integer> stringIntegerJavaPairRDD = sample.mapToPair(read -> new Tuple2<>(Utils.groupKey(read), 1)).reduceByKey((v1, v2) -> v1 + v2);
        System.out.println("total reads: " + toCache.count());
        System.out.println("keys: " + stringIntegerJavaPairRDD.count());
        Map<Integer, Integer> hist = stringIntegerJavaPairRDD.mapToPair(stringIntegerTuple2 -> new Tuple2<>(stringIntegerTuple2._2(), 1)).reduceByKey((v1, v2) -> v1 + v2).collectAsMap();

        SortedSet<Integer> sortedKeys = new TreeSet<>(hist.keySet());
        for(Integer i : sortedKeys) {
            System.out.println(i + "," + hist.get(i));
        }
    }

    static HashSet<String> getPathologicalKeys(JavaRDD<GATKRead> firstReads, double sampleRate, int pathologicalThreshold) {
        JavaRDD<GATKRead> sample = firstReads.sample(false, sampleRate);
        JavaPairRDD<String, Integer> numReadsByKey = sample.mapToPair(read -> new Tuple2<>(Utils.groupKey(read), 1)).reduceByKey((v1, v2) -> v1 + v2);
        Map<String, Integer> pathologicalKeys = numReadsByKey.filter(v1 -> v1._2() > pathologicalThreshold).collectAsMap();

        for (Map.Entry<String, Integer> e : pathologicalKeys.entrySet()) {
            System.out.println("key, count: " + e.getKey() + "-" + e.getValue());
        }
        return Sets.newHashSet(pathologicalKeys.keySet());
    }

    static Map<String, Tuple2<Iterable<GATKRead>, Integer>> getPathologicalData(
            JavaRDD<GATKRead> readsRDD, HashSet<String> pathologicalKeys) {
        JavaRDD<GATKRead> pathologicalReads = readsRDD.filter(v1 -> pathologicalKeys.contains(Utils.groupKey(v1)));
        pathologicalReads.cache();

        // Get the count and the representative read(s) for each key
        JavaRDD<GATKRead> dupes = pathologicalReads.filter(GATKRead::isDuplicate);
        if (pathologicalReads.isEmpty()) {
            return Maps.newHashMap();
        }
        HashMap<String, Integer> pathologicalDupCounts = dupes.map(v1 -> {
            HashMap<String, Integer> m = Maps.newHashMap();
            m.put(Utils.groupKey(v1), 1);
            return m;
        }).reduce((v1, v2) -> {
            for (Map.Entry<String, Integer> e : v2.entrySet()) {
                v1.put(e.getKey(), e.getValue());
            }
            return v1;
        });
        Map<String, Iterable<GATKRead>> representativeReads = pathologicalReads.filter(v1 -> !v1.isDuplicate()).mapToPair(read -> new Tuple2<>(Utils.groupKey(read), read)).groupByKey().collectAsMap();

        if (Sets.difference(pathologicalDupCounts.keySet(), representativeReads.keySet()).size() != 0) {
            throw new GATKException("Pathogogical sets should have the same keys");
        }
        Map<String, Tuple2<Iterable<GATKRead>, Integer>> pathologicalData = Maps.newHashMap();
        for (Map.Entry<String, Integer> e : pathologicalDupCounts.entrySet()) {
            pathologicalData.put(e.getKey(), new Tuple2<>(representativeReads.get(e.getKey()), e.getValue()));
        }
        pathologicalReads.unpersist();
        return pathologicalData;
    }


    public static Map<String, Integer> getPathologicalTagCounts(Map<String, Tuple2<Iterable<GATKRead>, Integer>> firstPathologicalResults,
                                                                Map<String, Tuple2<Iterable<GATKRead>, Integer>> secondPathologicalResults) {
        if (Sets.difference(firstPathologicalResults.keySet(), secondPathologicalResults.keySet()).size() != 0) {
            throw new GATKException("Pathological results should have the same keys");
        }

        List<String> pathologicalTags = Collections.emptyList();
        for (Map.Entry<String, Tuple2<Iterable<GATKRead>, Integer>> e : firstPathologicalResults.entrySet()) {
            Tuple2<Iterable<GATKRead>, Integer> firstT2 = e.getValue();
            Tuple2<Iterable<GATKRead>, Integer> secondT2 = secondPathologicalResults.get(e.getKey());
            if (firstT2._2().equals(secondT2._2())) {
                pathologicalTags.add("DIFF_NUM_DUPES");
            } else {
                Set<GATKRead> firstRepresentative = Sets.newHashSet(firstT2._1());
                Set<GATKRead> secondRepresentative = Sets.newHashSet(secondT2._1());
                if (Sets.difference(firstRepresentative, secondRepresentative).size() != 0) {
                    pathologicalTags.add("DUPES_MISMATCH");
                } else {
                    pathologicalTags.add("EQUAL");
                }
            }
        }
        Map<String, Integer> pathologicalTagCounts = Maps.newHashMap();

        pathologicalTags.forEach(s -> {
            Integer count = pathologicalTagCounts.getOrDefault(s, 0);
            pathologicalTagCounts.put(s, count+1);
        });
        return pathologicalTagCounts;
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

