package org.broadinstitute.hellbender.tools.spark.validation;


import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import htsjdk.samtools.SAMFileHeader;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.api.java.function.Function;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.TestSparkProgramGroup;
import org.broadinstitute.hellbender.engine.TraversalParameters;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.engine.spark.datasources.ReadsSparkSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadCoordinateComparator;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.read.markduplicates.ReadsKey;
import scala.Tuple2;

import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * The goal of this program is to load two potentially identical BAMs and determine if the BAMs contain the same
 * reads marked as duplicates. It fails if the BAMs don't have the same number of reads or total duplicates.
 */
@CommandLineProgramProperties(summary = "Compares two BAMs for duplicates", oneLineSummary = "Compares two BAMs for duplicates",
        programGroup = TestSparkProgramGroup.class)
@BetaFeature
public final class CompareDuplicatesSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;

    @Override
    public boolean requiresReads() { return true; }

    @Argument(doc="the second BAM", shortName = "I2", fullName = "input2", optional = false)
    protected String input2;

    @Argument(doc="print summary", shortName = "ps", fullName = "printSummary", optional = true)
    protected boolean printSummary = true;

    @Argument(doc="throw error on diff", shortName = "cd", fullName = "throwOnDiff", optional = true)
    protected boolean throwOnDiff = false;

    @Override
    protected void runTool(final JavaSparkContext ctx) {

        JavaRDD<GATKRead> firstReads = filteredReads(getReads(), readArguments.getReadFilesNames().get(0));

        ReadsSparkSource readsSource2 = new ReadsSparkSource(ctx, readArguments.getReadValidationStringency());
        TraversalParameters traversalParameters;
        if ( hasIntervals() ) {
            traversalParameters = intervalArgumentCollection.getTraversalParameters(getHeaderForReads().getSequenceDictionary());
        } else {
            traversalParameters = null;
        }
        JavaRDD<GATKRead> secondReads =  filteredReads(readsSource2.getParallelReads(input2, null, traversalParameters, bamPartitionSplitSize), input2);

        // Start by verifying that we have same number of reads and duplicates in each BAM.
        long firstBamSize = firstReads.count();

        long secondBamSize = secondReads.count();

        if (firstBamSize != secondBamSize) {
            throw new UserException("input bams have different numbers of mapped reads: "
                    + firstBamSize + "," + secondBamSize);
        }
        System.out.println("processing bams with " + firstBamSize + " mapped reads");


        long firstDupesCount = firstReads.filter(GATKRead::isDuplicate).count();
        long secondDupesCount = secondReads.filter(GATKRead::isDuplicate).count();
        if (firstDupesCount != secondDupesCount) {
            System.out.println("BAMs have different number of total duplicates: " + firstDupesCount + "," + secondDupesCount);
        }
        System.out.println("first and second: " + firstDupesCount + "," + secondDupesCount);

        Broadcast<SAMFileHeader> bHeader = ctx.broadcast(getHeaderForReads());
        // Group the reads of each BAM by MarkDuplicates key, then pair up the the reads for each BAM.
        JavaPairRDD<String, GATKRead> firstKeyed = firstReads.mapToPair(read -> new Tuple2<>(ReadsKey.keyForFragment(bHeader.getValue(), read), read));
        JavaPairRDD<String, GATKRead> secondKeyed = secondReads.mapToPair(read -> new Tuple2<>(ReadsKey.keyForFragment(bHeader.getValue(), read), read));
        JavaPairRDD<String, Tuple2<Iterable<GATKRead>, Iterable<GATKRead>>> cogroup = firstKeyed.cogroup(secondKeyed, getRecommendedNumReducers());


        // Produces an RDD of MatchTypes, e.g., EQUAL, DIFFERENT_REPRESENTATIVE_READ, etc. per MarkDuplicates key,
        // which is approximately start position x strand.
        JavaRDD<MatchType> tagged = cogroup.map(v1 -> {
            SAMFileHeader header = bHeader.getValue();

            Iterable<GATKRead> iFirstReads = v1._2()._1();
            Iterable<GATKRead> iSecondReads = v1._2()._2();

            return getDupes(iFirstReads, iSecondReads, header);
        });

        // TODO: We should also produce examples of reads that don't match to make debugging easier (#1263).
        Map<MatchType, Integer> tagCountMap = tagged.mapToPair(v1 ->
                new Tuple2<>(v1, 1)).reduceByKey((v1, v2) -> v1 + v2).collectAsMap();

        if (tagCountMap.get(MatchType.SIZE_UNEQUAL) != null) {
            throw new UserException("The number of reads by the MarkDuplicates key were unequal, indicating that the BAMs are not the same");
        }
        if (tagCountMap.get(MatchType.READ_MISMATCH) != null) {
            throw new UserException("The reads grouped by the MarkDuplicates key were not the same, indicating that the BAMs are not the same");
        }
        if (printSummary) {
            MatchType[] values = MatchType.values();
            Set<MatchType> matchTypes = Sets.newLinkedHashSet(Sets.newHashSet(values));
            System.out.println("##############################");
            matchTypes.forEach(s -> System.out.println(s + ": " + tagCountMap.getOrDefault(s, 0)));
        }

        if (throwOnDiff) {
            for (MatchType s : MatchType.values()) {
                if (s != MatchType.EQUAL) {
                    if (tagCountMap.get(s) != null)
                    throw new UserException("found difference between the two BAMs: " + s + " with count " + tagCountMap.get(s));
                }
            }
        }
    }

    /**
     * getDupes returns the metadata about how well the two sets of reads match.
     * @param f reads from the first bam
     * @param s reads from the second bam
     * @param header header (should be the same for both)
     * @return the type of the match, EQUAL, DIFFERENT_REPRESENTATIVE_READ, etc.
     */
    static MatchType getDupes(Iterable<GATKRead> f, Iterable<GATKRead> s, SAMFileHeader header) {
        List<GATKRead> first = Lists.newArrayList(f);
        List<GATKRead> second = Lists.newArrayList(s);
        if (first.size() != second.size()) {
            return MatchType.SIZE_UNEQUAL;
        }
        int size = first.size();
        first.sort(new ReadCoordinateComparator(header));
        second.sort(new ReadCoordinateComparator(header));

        Set<GATKRead> firstDupes = Sets.newLinkedHashSet();
        Set<GATKRead> secondDupes = Sets.newLinkedHashSet();
        for (int i = 0; i < size; ++i) {
            GATKRead firstRead = first.get(i);
            GATKRead secondRead = second.get(i);
            if (!firstRead.getName().equals(secondRead.getName())) {
                return MatchType.READ_MISMATCH;
            }
            if (firstRead.isDuplicate()) {
                firstDupes.add(firstRead);
            }
            if (secondRead.isDuplicate()) {
                secondDupes.add(secondRead);
            }
        }
        if (firstDupes.size() != secondDupes.size()) {
            return MatchType.DIFF_NUM_DUPES;
        }
        if (!firstDupes.equals(secondDupes)) {
            return  MatchType.DIFFERENT_REPRESENTATIVE_READ;
        }
        return  MatchType.EQUAL;
    }


    static JavaRDD<GATKRead> filteredReads(JavaRDD<GATKRead> initialReads, String fileName) {
        // We only need to compare duplicates that are "primary" (i.g., primary mapped read).
        return initialReads.map((Function<GATKRead, GATKRead>) v1 -> {
            v1.clearAttributes();
            return v1;
        }).filter(v1 -> {
            if (ReadUtils.isNonPrimary(v1) && v1.isDuplicate()) {
                throw new GATKException("found a non-primary read marked as a duplicate in the bam: "
                        + fileName);
            }
            return !ReadUtils.isNonPrimary(v1);
        });
    }

    enum MatchType {
        EQUAL,          // The same number of reads per key, same number of duplicates, and same representative read.
        SIZE_UNEQUAL,   // Difference number of reads at this key (this means the BAMs are different).
        READ_MISMATCH,  // The reads at this key are different for the two BAMs (which means the BAMs are different).
        DIFF_NUM_DUPES, // The reads are the same from each BAM, but there are different number of duplicates marked.
        DIFFERENT_REPRESENTATIVE_READ,  // The reads are the same and the same number of duplicates, but with a
        // different representative read (usually this is from a tie-breaker).
    }

}


