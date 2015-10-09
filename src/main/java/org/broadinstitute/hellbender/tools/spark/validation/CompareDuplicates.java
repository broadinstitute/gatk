package org.broadinstitute.hellbender.tools.spark.validation;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import htsjdk.samtools.SAMFileHeader;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.api.java.function.Function;
import org.apache.spark.api.java.function.PairFunction;
import org.apache.spark.api.java.function.VoidFunction;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.SparkProgramGroup;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.engine.spark.datasources.ReadsSparkSource;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.GoogleGenomicsReadToGATKReadAdapter;
import org.broadinstitute.hellbender.utils.read.ReadCoordinateComparator;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;
import scala.Tuple2;
import scala.Tuple3;

import java.util.HashSet;
import java.util.List;
import java.util.Map;

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

        JavaPairRDD<String, Iterable<GATKRead>> firstKeyed = firstReads.groupBy(v1 -> v1.getContig() + "," + v1.getStart());
        //JavaPairRDD<Long, Iterable<GATKRead>> secondKeyed = secondReads.groupBy(v1 -> (long) v1.getStart());
        JavaPairRDD<String, Iterable<GATKRead>> secondKeyed = secondReads.groupBy(v1 -> v1.getContig() + "," + v1.getStart());

        JavaPairRDD<String, Tuple2<Iterable<GATKRead>, Iterable<GATKRead>>> join = firstKeyed.join(secondKeyed);
        JavaPairRDD<String, Tuple3<String, Iterable<GATKRead>, Iterable<GATKRead>>> tagged = join.mapToPair(v1 -> {
            SAMFileHeader header = broadcastHeader.getValue();
            Iterable<GATKRead> firstReads1 = v1._2()._1();
            Iterable<GATKRead> secondReads1 = v1._2()._2();
            String str = new Dup().call(new Tuple3<>(firstReads1, secondReads1, header));
            return new Tuple2<>(v1._1(), new Tuple3<>(str, firstReads1, secondReads1));
        });
        JavaPairRDD<String, Tuple3<String, Iterable<GATKRead>, Iterable<GATKRead>>> equal =
                tagged.filter(v1 -> v1._2()._1().equals("EQUAL"));
        JavaPairRDD<String, Tuple3<String, Iterable<GATKRead>, Iterable<GATKRead>>> unequal =
                tagged.filter(v1 -> !v1._2()._1().equals("EQUAL"));

        Map<String, Tuple3<String, Iterable<GATKRead>, Iterable<GATKRead>>> diffNumDupesMap = unequal.filter(v1 -> v1._2()._1().equals("DIFF_NUM_DUPES")).collectAsMap();
        long diffNumDupes = diffNumDupesMap.size();
        long dupesMismatch = unequal.filter(v1 -> v1._2()._1().equals("DUPES_MISMATCH")).count();

        JavaRDD<GATKRead> firstDupes = firstReads.filter(GATKRead::isDuplicate);
        JavaRDD<GATKRead> secondDupes = secondReads.filter(GATKRead::isDuplicate);

        JavaRDD<GATKRead> firstDistinct = firstDupes.distinct();
        JavaRDD<GATKRead> secondDistinct = secondDupes.distinct();

        long allDupes = firstDistinct.union(secondDistinct).distinct().count();
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
        System.out.println("first and second, total distinct: " + firstDistinct.count() + "," +
                secondDistinct.count() + "," + allDupes);
        System.out.println("no diff,diffs,diffNumDupes,dupesMismatch: " + equal.count() + "," + unequal.count() + "," +
                diffNumDupes + "," + dupesMismatch);


        for (Map.Entry<String, Tuple3<String, Iterable<GATKRead>, Iterable<GATKRead>>> entry : diffNumDupesMap.entrySet()) {
            System.out.println("**************************************");
            System.out.println("*************** " + entry.getKey() + "*************************");
            System.out.println("**************************************");
            List<GATKRead> iFirst = Lists.newArrayList(entry.getValue()._2());
            List<GATKRead> iSecond = Lists.newArrayList(entry.getValue()._3());
            System.out.println("first");
            for (GATKRead r : iFirst) {
                System.out.println(((SAMRecordToGATKReadAdapter) r).getSamRecord().getSAMString());
            }
            System.out.println("second");
            for (GATKRead r : iSecond) {
                System.out.println(((SAMRecordToGATKReadAdapter) r).getSamRecord().getSAMString());
            }
        }
        

        /*
        unequal.foreach(v1 -> System.out.println(v1._2()._1() + "," +
                Lists.newArrayList(v1._2()._2()).size() + "," +
                Lists.newArrayList(v1._2()._3()).size()));
        */
    }
}

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
}

