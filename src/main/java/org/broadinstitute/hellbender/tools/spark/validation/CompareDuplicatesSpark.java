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
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.cmdline.programgroups.DiagnosticsAndQCProgramGroup;
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
 * Determine if two potentially identical BAMs have the same duplicate reads. This tool is useful for checking if two
 * BAMs that seem identical have the same reads marked as duplicates.
 *
 * <p>It fails if at least one of the following conditions is true of the two BAM files:</p>
 *
 * <ul>
 *     <li>Different number of primary mapped reads</li>
 *     <li>Different number of duplicate reads (as indicated by the SAM record flag)</li>
 *     <li>Different reads mapped to some position in the reference</li>
 * </ul>
 *
 * <p>The tool gathers the mapped reads together into groups that belong to the same library and map to the same
 * position and strand in the reference. If the tool does not fail, then it reports the number of these groups with
 * the following properties:</p>
 *
 * <ul>
 *     <li>SIZE_UNEQUAL: different number of reads</li>
 *     <li>EQUAL: same reads and same duplicates</li>
 *     <li>READ_MISMATCH: reads with different names</li>
 *     <li>DIFFERENT_REPRESENTATIVE_READ: same reads and number of duplicates, but one or more duplicate reads are different</li>
 *     <li>DIFF_NUM_DUPES: same reads but different number of duplicates</li>
 * </ul>
 *
 * <h3>Input</h3>
 *
 * <ul>
 *     <li>Two BAM files</li>
 * </ul>
 *
 * <h3>Output</h3>
 *
 * <p>Results are printed to the screen</p>
 *
 * <h3>Usage example</h3>
 * <pre>
 * gatk CompareDuplicatesSpark \
 *     -I input_reads_1.bam \
 *     -I2 input_reads_2.bam
 * </pre>
 *
 * <p>This tool can be run without explicitly specifying Spark options. That is to say, the given example command
 * without Spark options will run locally. See
 * <a href ="https://software.broadinstitute.org/gatk/documentation/article?id=10060">Tutorial#10060</a> for an example
 * of how to set up and run a Spark tool on a cloud Spark cluster.</p>
 */
@DocumentedFeature
@CommandLineProgramProperties(summary = "Determine if two potentially identical BAMs have the same duplicate reads. " +
        "This tool is useful for checking if two BAMs that seem identical have the same reads marked as duplicates.",
        oneLineSummary = "Determine if two potentially identical BAMs have the same duplicate reads",
        programGroup = DiagnosticsAndQCProgramGroup.class)
@BetaFeature
public final class CompareDuplicatesSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;

    public static final String INPUT_2_LONG_NAME = "input2";
    public static final String INPUT_2_SHORT_NAME = "I2";
    public static final String PRINT_SUMMARY_LONG_NAME = "print-summary";
    public static final String THROW_ON_DIFF_LONG_NAME = "throw-on-diff";

    @Override
    public boolean requiresReads() { return true; }

    @Argument(doc="The second BAM", shortName = INPUT_2_SHORT_NAME, fullName = INPUT_2_LONG_NAME, optional = false)
    protected String input2;

    @Argument(doc="Print a summary", fullName = PRINT_SUMMARY_LONG_NAME, optional = true)
    protected boolean printSummary = true;

    @Argument(doc="Throw error if any differences were found", fullName = THROW_ON_DIFF_LONG_NAME, optional = true)
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
        //TODO these are probably horribly broken forever and they should be checked given the new key scheme
        JavaPairRDD<Integer, GATKRead> firstKeyed = firstReads.mapToPair(read -> new Tuple2<>(ReadsKey.hashKeyForFragment(
                ReadUtils.getStrandedUnclippedStart(
                                                                                                                                  read),
                read.isReverseStrand(),
                ReadUtils.getReferenceIndex(
                                                                                                                                  read,
                                                                                                                                  bHeader.getValue()),
                ReadUtils.getLibrary(read, bHeader.getValue())), read));
        JavaPairRDD<Integer, GATKRead> secondKeyed = secondReads.mapToPair(read -> new Tuple2<>(ReadsKey.hashKeyForFragment(
                ReadUtils.getStrandedUnclippedStart(
                                                                                                                                    read),
                read.isReverseStrand(),
                ReadUtils.getReferenceIndex(
                                                                                                                                    read,
                                                                                                                                    bHeader.getValue()),
                ReadUtils.getLibrary(read, bHeader.getValue())), read));
        JavaPairRDD<Integer, Tuple2<Iterable<GATKRead>, Iterable<GATKRead>>> cogroup = firstKeyed.cogroup(secondKeyed, getRecommendedNumReducers());


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


