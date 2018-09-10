package org.broadinstitute.hellbender.tools.walkers.markduplicates;

import htsjdk.samtools.*;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import org.apache.commons.io.FileUtils;
import org.apache.spark.SparkException;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.ReadsDataSource;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.BaseTest;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.broadinstitute.hellbender.testutils.SamAssertionUtils;
import org.broadinstitute.hellbender.testutils.testers.MarkDuplicatesSparkTester;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.markduplicates.GATKDuplicationMetrics;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.sam.util.PhysicalLocationInt;
import picard.sam.util.ReadNameParser;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

/**
 * This class defines the individual test cases to run.  The actual running of the test is done
 * by AbstractMarkDuplicatesCommandLineProgramTester or children thereof (see getTester).
 */
public abstract class AbstractMarkDuplicatesCommandLineProgramTest extends CommandLineProgramTest {

    public static final File TEST_DATA_DIR = new File(getTestDataDir(), "walkers/MarkDuplicatesGATK/");

    protected abstract MarkDuplicatesSparkTester getTester();

    protected abstract CommandLineProgram getCommandLineProgramInstance();

    protected boolean markSecondaryAndSupplementaryRecordsLikeTheCanonical() { return false; }

    // ELIGIBLE_BASE_QUALITY is the minimum quality considered by the dataflow version
    // for marking a specific pair as best. We use it to ensure that the expected
    // pair wins the tie so that the tie is not broken randomly.
    protected static final int DEFAULT_BASE_QUALITY = 10;
    protected static final int ELIGIBLE_BASE_QUALITY = 15;

    @Test
    public void testSingleUnmappedFragment() {
        final MarkDuplicatesSparkTester tester = getTester();
        tester.addUnmappedFragment(-1, DEFAULT_BASE_QUALITY);
        tester.runTest();
    }

    @Test
    public void testTwoUnmappedFragments() {
        final MarkDuplicatesSparkTester tester = getTester();
        tester.addUnmappedFragment(-1, DEFAULT_BASE_QUALITY);
        tester.addUnmappedFragment(-1, DEFAULT_BASE_QUALITY);
        tester.runTest();
    }

    @Test
    public void testSingleUnmappedPair() {
        final MarkDuplicatesSparkTester tester = getTester();
        tester.addUnmappedPair(-1, DEFAULT_BASE_QUALITY);
        tester.runTest();
    }

    @Test
    public void testSingleMappedFragment() {
        final MarkDuplicatesSparkTester tester = getTester();
        tester.addMappedFragment(1, 1, false, DEFAULT_BASE_QUALITY);
        tester.runTest();
    }

    @Test
    public void testTwoMappedFragments() {
        final MarkDuplicatesSparkTester tester = getTester();
        tester.addMappedFragment(0, 1, false, ELIGIBLE_BASE_QUALITY);
        tester.addMappedFragment(0, 1, true, DEFAULT_BASE_QUALITY); // duplicate!!!
        tester.runTest();
    }

    @Test
    public void testSingleMappedPair() {
        final MarkDuplicatesSparkTester tester = getTester();
        tester.addMappedPair(1, 1, 100, false, false, DEFAULT_BASE_QUALITY);
        tester.runTest();
    }

    @Test
    public void testSingleMappedFragmentAndSingleMappedPair() {
        final MarkDuplicatesSparkTester tester = getTester();
        tester.addMappedFragment(1, 1, true, DEFAULT_BASE_QUALITY); // duplicate!!!
        tester.addMappedPair(1, 1, 100, false, false, ELIGIBLE_BASE_QUALITY);
        tester.runTest();
    }

    @Test
    public void testTwoMappedPairs() {
        final MarkDuplicatesSparkTester tester = getTester();
        tester.addMappedPair(1, 1, 100, false, false, ELIGIBLE_BASE_QUALITY);
        tester.addMappedPair(1, 1, 100, true, true, DEFAULT_BASE_QUALITY); // duplicate!!!
        tester.runTest();
    }

    @Test
    public void testThreeMappedPairs() {
        final MarkDuplicatesSparkTester tester = getTester();
        tester.addMappedPair(1, 1, 100, false, false, ELIGIBLE_BASE_QUALITY);
        tester.addMappedPair(1, 1, 100, true, true, DEFAULT_BASE_QUALITY); // duplicate!!!
        tester.addMappedPair(1, 1, 100, true, true, DEFAULT_BASE_QUALITY); // duplicate!!!
        tester.runTest();
    }

    @Test
    public void testSingleMappedFragmentAndTwoMappedPairs() {
        final MarkDuplicatesSparkTester tester = getTester();
        tester.addMappedFragment(1, 1, true, DEFAULT_BASE_QUALITY); // duplicate!!!
        tester.addMappedPair(1, 1, 100, false, false, ELIGIBLE_BASE_QUALITY);
        tester.addMappedPair(1, 1, 100, true, true, DEFAULT_BASE_QUALITY); // duplicate!!!
        tester.runTest();
    }

    @Test
    public void testTwoMappedPairsAndTerminalUnmappedFragment() {
        final MarkDuplicatesSparkTester tester = getTester();
        tester.addMappedPair(1, 1, 100, false, false, ELIGIBLE_BASE_QUALITY);
        tester.addMappedPair(1, 1, 100, true, true, DEFAULT_BASE_QUALITY); // duplicate!!!
        tester.addUnmappedFragment(-1, DEFAULT_BASE_QUALITY); // unmapped fragment at end of file
        tester.runTest();
    }

    @Test
    public void testTwoMappedPairsAndTerminalUnmappedPair() {
        final MarkDuplicatesSparkTester tester = getTester();
        tester.addMappedPair(1, 1, 100, false, false, ELIGIBLE_BASE_QUALITY);
        tester.addMappedPair(1, 1, 100, true, true, DEFAULT_BASE_QUALITY); // duplicate!!!
        tester.addUnmappedPair(-1, DEFAULT_BASE_QUALITY); // unmapped pair at end of file
        tester.runTest();
    }

    @Test
    public void testOpticalDuplicateFinding() {
        final MarkDuplicatesSparkTester tester = getTester();

        // explicitly creating 1 expected optical duplicate pair
        tester.setExpectedOpticalDuplicate(1);

        // pass in the read names manually, in order to control duplicates vs optical duplicates
        tester.addMatePair("READ0:1:1:1:1", 1, 1, 100, false, false, false, false, "50M", "50M", false, true, false,
                           false, false, ELIGIBLE_BASE_QUALITY); // non-duplicate mapped pair to start
        tester.addMatePair("READ1:1:1:1:300", 1, 1, 100, false, false, true, true, "50M", "50M", false, true, false,
                           false, false, DEFAULT_BASE_QUALITY); // duplicate pair, NOT optical duplicate (delta-Y > 100)
        tester.addMatePair("READ2:1:1:1:50", 1, 1, 100, false, false, true, true, "50M", "50M", false, true, false,
                           false, false, DEFAULT_BASE_QUALITY); // duplicate pair, expected optical duplicate (delta-X and delta-Y < 100)
        tester.runTest();
    }

    @Test
    public void testOpticalDuplicatesDifferentReadGroups() {
        final MarkDuplicatesSparkTester tester = getTester();
        tester.setExpectedOpticalDuplicate(0);
        tester.addMatePair("RUNID:7:1203:2886:82292",  19, 19, 485253, 485253, false, false, true, true, "42M59S", "59S42M", true, false, false, false, false, DEFAULT_BASE_QUALITY, "H0164.2");  // duplicate
        tester.addMatePair("RUNID:7:1203:2886:82242", 19, 19, 485253, 485253, false, false, false, false, "42M59S", "59S42M", true, false, false, false, false, DEFAULT_BASE_QUALITY, "H0164.1");
        SAMFileHeader header = tester.getHeader();
        SAMReadGroupRecord readGroup1 = new SAMReadGroupRecord("H0164.2");
        SAMReadGroupRecord readGroup2 = new SAMReadGroupRecord("H0164.1");
        readGroup1.setSample("test");
        readGroup2.setSample("test");
        header.addReadGroup(readGroup1);
        header.addReadGroup(readGroup2);
        tester.runTest();
    }

    @Test
    public void testOpticalDuplicatesTheSameReadGroup() {
        final MarkDuplicatesSparkTester tester = getTester();
        tester.setExpectedOpticalDuplicate(1);
        tester.addMatePair("RUNID:7:1203:2886:82292",  19, 19, 485253, 485253, false, false, true, true, "42M59S", "59S42M", true, false, false, false, false, DEFAULT_BASE_QUALITY, "H0164.2");  // duplicate (tie broken by readname)
        tester.addMatePair("RUNID:7:1203:2886:82242", 19, 19, 485253, 485253, false, false, false, false, "42M59S", "59S42M", true, false, false, false, false, DEFAULT_BASE_QUALITY, "H0164.2");

        SAMFileHeader header = tester.getHeader();
        SAMReadGroupRecord readGroup1 = new SAMReadGroupRecord("H0164.2");
        SAMReadGroupRecord readGroup2 = new SAMReadGroupRecord("H0164.1");
        readGroup1.setSample("test");
        readGroup2.setSample("test");
        header.addReadGroup(readGroup1);
        header.addReadGroup(readGroup2);
        tester.runTest();
    }

    @Test
    public void testOpticalDuplicatesAndPCRDuplicatesOrientationDifference() {
        final MarkDuplicatesSparkTester tester = getTester();
        tester.setExpectedOpticalDuplicate(0);
        tester.addMatePair("RUNID:7:1203:2886:16838",  19, 19, 485253, 486253, false, false, true, true, "101M", "101M", true, false, false, false, false, DEFAULT_BASE_QUALITY, "1");  // duplicate
        tester.addMatePair("RUNID:7:1203:2886:16834", 19, 19, 486253, 485253, false, false, false, false, "101M", "101M", false, true, false, false, false, DEFAULT_BASE_QUALITY, "1");
        // Even though these reads are in a duplicate group together, we don't want to mark them as Optical Duplicates because their orientation is flipped (which doesn't matter for PCR duplicates)
        tester.runTest();
    }

    @Test
    // To match picard, when two reads are grouped into a duplicate group and have the same unclipped start position we will unify their orientations
    // to FR if they point in opposite directions. This test asserts that two pairs of reads that ultimately all map to 93470499-93470499 with clipping
    // and orientations will ultimately be marked and mapped together even though they have opposite orientations (one is FR while the other is RF).
    public void testReadSameStartPositionOrientationOverride() {
        final MarkDuplicatesSparkTester tester = getTester();
        tester.setExpectedOpticalDuplicate(0);
        tester.addMatePair("HJYFJCCXX160204:8:1124:31659:21526",  19, 19, 93470380, 93470499, false, false, false, false, "31S111M9S", "64M87S", true, false, false, false, false, DEFAULT_BASE_QUALITY, "1");  // after clipping these both start on 93470499
        tester.addMatePair("HJYFJCCXX160204:8:1124:31659:21527",  19, 19, 93470499, 93470380, false, false, true, true, "64M87S", "31S111M9S", false, true, false, false, false, DEFAULT_BASE_QUALITY, "1");  // after clipping these both start on 93470499
        tester.runTest();
    }

    @Test
    // This asserts that we are flipping reads in the Pair object based on both start position and contig index, this does not
    // make a difference unless the start position is the same across two contigs, so we assert it is handled properly
    public void testReadsHaveSameStartPositionButDifferentChromosomeNonEquivalence() {
        final MarkDuplicatesSparkTester tester = getTester();
        tester.setExpectedOpticalDuplicate(0);
        tester.addMatePair("HJYFJCCXX160204:8:1124:31659:21526",  19, 20, 93470380, 93470499, false, false, false, false, "31S111M9S", "64M87S", true, false, false, false, false, DEFAULT_BASE_QUALITY, "1");  // after clipping these both start on 93470499
        tester.addMatePair("HJYFJCCXX160204:8:1124:31659:21527",  20, 19, 93470499, 93470380, false, false, true, true, "64M87S", "31S111M9S", false, true, false, false, false, DEFAULT_BASE_QUALITY, "1");  // after clipping these both start on 93470499
        tester.runTest();
    }

    @Test
    public void testOpticalDuplicateClusterSamePositionNoOpticalDuplicatesWithinPixelDistance() {
        final MarkDuplicatesSparkTester tester = getTester();
        tester.getSamRecordSetBuilder().setReadLength(101);
        tester.setExpectedOpticalDuplicate(0);
        tester.addMatePair("RUNID:7:1203:2886:16834", 1, 485253, 485253, false, false, true, true, "42M59S", "59S42M", false, true, false, false, false, DEFAULT_BASE_QUALITY); // duplicate
        tester.addMatePair("RUNID:7:1203:2884:16835", 1, 485253, 485253, false, false, false, false, "59S42M", "42M59S", true, false, false, false, false, ELIGIBLE_BASE_QUALITY);
        tester.runTest();
    }

    @Test
    public void testOpticalDuplicateClusterSamePositionOneOpticalDuplicatesWithinPixelDistance() {
        final MarkDuplicatesSparkTester tester = getTester();
        tester.getSamRecordSetBuilder().setReadLength(45);
        tester.setExpectedOpticalDuplicate(1);
        tester.addMatePair("RUNID:7:1203:2886:16834", 1, 485253, 485253, false, false, true, true, "45M", "45M", false, true, false, false, false, DEFAULT_BASE_QUALITY); // duplicate
        tester.addMatePair("RUNID:7:1203:2884:16835", 1, 485253, 485253, false, false, false, false, "45M", "45M", false, true, false, false, false, ELIGIBLE_BASE_QUALITY);
        tester.runTest();
    }

    @Test
    public void testOpticalDuplicateClusterOneEndSamePositionOneCluster() {
        final MarkDuplicatesSparkTester tester = getTester();
        tester.getSamRecordSetBuilder().setReadLength(101);
        tester.setExpectedOpticalDuplicate(1);
        tester.addMatePair("RUNID:7:2205:17939:39728", 1, 485328, 485312, false, false, false, false, "55M46S", "30S71M", false, true, false, false, false, ELIGIBLE_BASE_QUALITY);
        tester.addMatePair("RUNID:7:2205:17949:39745", 1, 485328, 485328, false, false, true, true, "55M46S", "46S55M", false, true, false, false, false, DEFAULT_BASE_QUALITY); // duplicate
        tester.runTest();
    }

    @Test
    public void testTwoMappedPairsAndMappedSecondaryFragment() {
        final MarkDuplicatesSparkTester tester = getTester();
        tester.addMappedPair(1, 1, 100, false, false, ELIGIBLE_BASE_QUALITY);
        tester.addMappedPair(1, 1, 100, true, true, DEFAULT_BASE_QUALITY); // duplicate!!!
        tester.addMappedFragment(1, 200, false, DEFAULT_BASE_QUALITY, true); // mapped non-primary fragment
        tester.runTest();
    }

    @Test
    public void testMappedFragmentAndMappedPairFirstOfPairNonPrimary() {
        final MarkDuplicatesSparkTester tester = getTester();
        tester.getSamRecordSetBuilder().setReadLength(76);
        tester.addMappedFragment(1, 1, false, ELIGIBLE_BASE_QUALITY);
        tester.addMatePair(1, 200, 0, false, true, false, false, "54M22S", null, false, false, true, true, false, DEFAULT_BASE_QUALITY);
        tester.runTest();
    }

    @Test
    public void testTwoMappedPairsMatesSoftClipped() {
        final MarkDuplicatesSparkTester tester = getTester();
        tester.getSamRecordSetBuilder().setReadLength(76);
        tester.addMappedPair(1, 10022, 10051, false, false, "76M", "8S68M", false, true, false, DEFAULT_BASE_QUALITY);
        tester.addMappedPair(1, 10022, 10063, false, false, "76M", "5S71M", false, true, false, DEFAULT_BASE_QUALITY);
        tester.runTest();
    }


    @Test
    public void testTwoMappedPairsWithSoftClipping() {
        final MarkDuplicatesSparkTester tester = getTester();
        tester.getSamRecordSetBuilder().setReadLength(76);
        //todo Picard says no duplicates, we say duplicate, what's going on?
        // NB: no duplicates
        // 5'1: 2, 5'2:46+73M=118
        // 5'1: 2, 5'2:51+68M=118
        tester.addMappedPair(1, 2, 46, false, false, "6S42M28S", "3S73M", false, ELIGIBLE_BASE_QUALITY);
        tester.addMappedPair(1, 2, 51, true, true, "6S42M28S", "8S68M", false, DEFAULT_BASE_QUALITY); // duplicate
        tester.runTest();
    }

    @Test
    public void testTwoMappedPairsWithSoftClippingFirstOfPairOnlyNoMateCigar() {
        final MarkDuplicatesSparkTester tester = getTester();
        tester.getSamRecordSetBuilder().setReadLength(76);
        tester.setNoMateCigars(true);
        // NB: no duplicates
        // 5'1: 2, 5'2:46+73M=118
        // 5'1: 2, 5'2:51+68M=118
        tester.addMappedPair(1, 12, 46, false, false, "6S42M28S", null, true, DEFAULT_BASE_QUALITY); // only add the first one
        tester.addMappedPair(1, 12, 51, false, false, "6S42M28S", null, true, DEFAULT_BASE_QUALITY); // only add the first one
        tester.runTest();
    }

    @Test
    public void testTwoMappedPairsWithSoftClippingBoth() {
        final MarkDuplicatesSparkTester tester = getTester();
        tester.getSamRecordSetBuilder().setReadLength(76);
        // mapped reference length: 73 + 42 = 115
        tester.addMappedPair(1, 10046, 10002, true, true, "3S73M", "6S42M28S", true, false, false, DEFAULT_BASE_QUALITY); // duplicate
        // mapped reference length: 68 + 48 = 116
        tester.addMappedPair(1, 10051, 10002, false, false, "8S68M", "6S48M22S", true, false, false, ELIGIBLE_BASE_QUALITY);
        tester.runTest();
    }

    @Test
    public void testMatePairSecondUnmapped() {
        final MarkDuplicatesSparkTester tester = getTester();
        tester.getSamRecordSetBuilder().setReadLength(76);
        tester.addMatePair(1, 10049, 10049, false, true, false, false, "11M2I63M", null, false, false, false, false, false, DEFAULT_BASE_QUALITY);   // neither are duplicates
        tester.runTest();
    }

    @Test
    public void testMatePairFirstUnmapped() {
        final MarkDuplicatesSparkTester tester = getTester();
        tester.getSamRecordSetBuilder().setReadLength(76);
        tester.addMatePair(1, 10056, 10056, true, false, false, false, null, "54M22S", false, false, false, false, false, DEFAULT_BASE_QUALITY);    // neither are duplicates
        tester.runTest();
    }

    @Test
    public void testMappedFragmentAndMatePairSecondUnmapped() {
        final MarkDuplicatesSparkTester tester = getTester();
        tester.getSamRecordSetBuilder().setReadLength(76);
        tester.addMatePair(1, 10049, 10049, false, true, false, false, "11M2I63M", null, false, false, false, false, false, ELIGIBLE_BASE_QUALITY);
        // We set the length here separately because some of the MD variants use TOTAL_MAPPED_REFERENCE_LENGTH as their DUPLICATE_SCORING_STRATEGY.  If we kept the read length
        // as 76 like with the mate pair, it would mark the mapped read in the pair as a duplicate because it has less reference length due to the two insertions in the cigar
        tester.getSamRecordSetBuilder().setReadLength(36);
        tester.addMappedFragment(1, 10049, true, DEFAULT_BASE_QUALITY); // duplicate
        tester.runTest();
    }

    @Test
    public void testMappedFragmentAndMatePairFirstUnmapped() {
        final MarkDuplicatesSparkTester tester = getTester();
        tester.getSamRecordSetBuilder().setReadLength(76);
        tester.addMatePair(1, 10049, 10049, true, false, false, false, null, "11M2I63M", false, false, false, false, false, ELIGIBLE_BASE_QUALITY);
           // We set the length here separately because some of the MD variants use TOTAL_MAPPED_REFERENCE_LENGTH as their DUPLICATE_SCORING_STRATEGY.  If we kept the read length
        // as 76 like with the mate pair, it would mark the mapped read in the pair as a duplicate because it has less reference length due to the two insertions in the cigar
        tester.getSamRecordSetBuilder().setReadLength(36);
        tester.addMappedFragment(1, 10049, true, DEFAULT_BASE_QUALITY); // duplicate
        tester.runTest();
    }

    @Test
    public void testMappedPairAndMatePairSecondUnmapped() {
        final MarkDuplicatesSparkTester tester = getTester();
        tester.getSamRecordSetBuilder().setReadLength(76);
        tester.addMatePair(1, 10040, 10040, false, true, true, false, "76M", null, false, false, false, false, false, DEFAULT_BASE_QUALITY); // second a duplicate,
        // second end unmapped
        tester.addMappedPair(1, 10189, 10040, false, false, "41S35M", "65M11S", true, false, false, ELIGIBLE_BASE_QUALITY); // mapped OK
        tester.runTest();
    }

    @Test
    public void testMappedPairAndMatePairFirstUnmapped() {
        final MarkDuplicatesSparkTester tester = getTester();
        tester.getSamRecordSetBuilder().setReadLength(76);
        tester.addMatePair(1, 10040, 10040, true, false, false, true,  null, "76M", false, false, false, false, false, DEFAULT_BASE_QUALITY); // first a duplicate,
        // first end unmapped
        tester.addMappedPair(1, 10189, 10040, false, false, "41S35M", "65M11S", true, false, false, ELIGIBLE_BASE_QUALITY); // mapped OK
        tester.runTest();
    }

    // TODO: fails on MarkDuplicatesWithMateCigar
    @Test
    public void testMappedPairAndMatePairFirstOppositeStrandSecondUnmapped() {
        final MarkDuplicatesSparkTester tester = getTester();
        tester.getSamRecordSetBuilder().setReadLength(101);
        // first end mapped OK -, second end unmapped
        tester.addMatePair(1, 484071, 484071, false, true, false, false,  "66S35M", null, true, false, false, false, false, DEFAULT_BASE_QUALITY);
        // mapped OK +/-
        tester.addMappedPair(1, 484105, 484075, false, false, "35M66S", "30S71M", false, true, false, DEFAULT_BASE_QUALITY);
        tester.runTest();
    }

    @Test
    public void testSupplementaryReadUnmappedMate() {
        File output = createTempFile("supplementaryReadUnmappedMate", "bam");
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addOutput(output);
        args.addInput(getTestFile("supplementaryReadUnmappedmate.bam"));
        runCommandLine(args);

        try ( final ReadsDataSource outputReadsSource = new ReadsDataSource(output.toPath()) ) {
            final List<GATKRead> actualReads = new ArrayList<>();
            for ( final GATKRead read : outputReadsSource ) {
                Assert.assertFalse(read.isDuplicate());
                actualReads.add(read);
            }

            Assert.assertEquals(actualReads.size(), 3, "Wrong number of reads output");
        }
    }

    @Test
    public void testMappedPairAndMappedFragmentAndMatePairSecondUnmapped() {
        final MarkDuplicatesSparkTester tester = getTester();
        tester.getSamRecordSetBuilder().setReadLength(76);
        tester.addMatePair(1, 10040, 10040, false, true, true, false, "76M", null, false, false, false, false, false, DEFAULT_BASE_QUALITY); // first a duplicate,
        // second end unmapped
        tester.addMappedPair(1, 10189, 10040, false, false, "41S35M", "65M11S", true, false, false, ELIGIBLE_BASE_QUALITY); // mapped OK
        tester.addMappedFragment(1, 10040, true, DEFAULT_BASE_QUALITY); // duplicate
        tester.runTest();
    }

    @Test
    public void testMappedPairAndMappedFragmentAndMatePairFirstUnmapped() {
        final MarkDuplicatesSparkTester tester = getTester();
        tester.getSamRecordSetBuilder().setReadLength(76);
        tester.addMatePair(1, 10040, 10040, true, false, false, true, null, "76M", false, false, false, false, false, DEFAULT_BASE_QUALITY); // first a duplicate,
        // first end unmapped
        tester.addMappedPair(1, 10189, 10040, false, false, "41S35M", "65M11S", true, false, false, ELIGIBLE_BASE_QUALITY); // mapped OK
        tester.addMappedFragment(1, 10040, true, DEFAULT_BASE_QUALITY); // duplicate
        tester.runTest();
    }

    @Test
    public void testTwoMappedPairsWithOppositeOrientations() {
        final MarkDuplicatesSparkTester tester = getTester();
        tester.getSamRecordSetBuilder().setReadLength(76);
        tester.addMappedPair(1, 10182, 10038, true, true, "32S44M", "66M10S", true, false, false, DEFAULT_BASE_QUALITY); // -/+
        tester.addMappedPair(1, 10038, 10182, false, false, "70M6S", "32S44M", false, true, false, ELIGIBLE_BASE_QUALITY); // +/-, both are duplicates
        tester.runTest();
    }

    @Test
    public void testTwoMappedPairsWithOppositeOrientationsNumberTwo() {
        final MarkDuplicatesSparkTester tester = getTester();
        tester.getSamRecordSetBuilder().setReadLength(76);
        tester.addMappedPair(1, 10038, 10182, false, false, "70M6S", "32S44M", false, true, false, ELIGIBLE_BASE_QUALITY); // +/-, both are duplicates
        tester.addMappedPair(1, 10182, 10038, true, true, "32S44M", "66M10S", true, false, false, DEFAULT_BASE_QUALITY); // -/+
        tester.runTest();
    }

    @Test
    public void testThreeMappedPairsWithMatchingSecondMate() {
        final MarkDuplicatesSparkTester tester = getTester();
        tester.getSamRecordSetBuilder().setReadLength(76);
        // Read0 and Read2 are duplicates
        // 10181+41=10220, 10058
        tester.addMappedPair(1, 10181, 10058, false, false, "41S35M", "47M29S", true, false, false, ELIGIBLE_BASE_QUALITY); // -/+
        // 10181+37=10216, 10058
        tester.addMappedPair(1, 10181, 10058, true, true, "37S39M", "44M32S", true, false, false, DEFAULT_BASE_QUALITY); // -/+
        // 10180+36=10216, 10058
        tester.addMappedPair(1, 10180, 10058, false, false, "36S40M", "50M26S", true, false, false, ELIGIBLE_BASE_QUALITY); // -/+, both are duplicates
        tester.runTest();
    }

    @Test
    public void testMappedPairWithSamePosition() {
        final MarkDuplicatesSparkTester tester = getTester();
        tester.getSamRecordSetBuilder().setReadLength(76);
        tester.addMappedPair(1, 4914, 4914, false, false, "37M39S", "73M3S", false, false, false, DEFAULT_BASE_QUALITY); // +/+
        tester.runTest();
    }

    @Test
    public void testMappedPairWithSamePositionSameCigar() {
        final MarkDuplicatesSparkTester tester = getTester();
        tester.getSamRecordSetBuilder().setReadLength(76);
        tester.addMappedPair(1, 4914, 4914, false, false, "37M39S", "37M39S", false, false, false, DEFAULT_BASE_QUALITY); // +/+
        tester.runTest();
    }

    @Test
    public void testTwoMappedPairWithSamePosition() {
        final MarkDuplicatesSparkTester tester = getTester();
        tester.getSamRecordSetBuilder().setReadLength(76);
        tester.addMappedPair(0, 5604914, 5604914, false, false, "37M39S", "73M3S", false, false, false, ELIGIBLE_BASE_QUALITY); // +/+
        tester.addMappedPair(0, 5604914, 5604914, true, true, "37M39S", "73M3S", false, false, false, DEFAULT_BASE_QUALITY); // +/+
        tester.runTest();
    }

    @Test
    public void testTwoMappedPairWithSamePositionDifferentStrands() {
        final MarkDuplicatesSparkTester tester = getTester();
        tester.addMappedPair(0, 5604914, 5604914, false, false, "50M", "50M", true, false, false, ELIGIBLE_BASE_QUALITY); // +/-
        tester.addMappedPair(0, 5604914, 5604914, true, true, "50M", "50M", false, true, false, DEFAULT_BASE_QUALITY); // -/+
        tester.runTest();
    }

    @Test
    public void testTwoMappedPairWithSamePositionDifferentStrands2() {
        final MarkDuplicatesSparkTester tester = getTester();
        tester.addMappedPair(0, 5604914, 5604915, false, false, "50M", "50M", true, false, false, ELIGIBLE_BASE_QUALITY); // +/-
        tester.addMappedPair(0, 5604915, 5604914, true, true, "50M", "50M", false, true, false, DEFAULT_BASE_QUALITY); // -/+
        tester.runTest();
    }

    @Test
    public void testMappedPairWithFirstEndSamePositionAndOther() {
        final MarkDuplicatesSparkTester tester = getTester();
        tester.getSamRecordSetBuilder().setReadLength(76);
        tester.addMappedPair(0, 5604914, 5605914, false, false, "37M39S", "73M3S", false, false, false, DEFAULT_BASE_QUALITY); // +/+
        tester.addMappedPair(0, 5604914, 5604914, false, false, "37M39S", "73M3S", false, false, false, DEFAULT_BASE_QUALITY); // +/+
        tester.runTest();
    }

    @Test
    public void testTwoGroupsOnDifferentChromosomesOfTwoFragments() {
        final MarkDuplicatesSparkTester tester = getTester();
        tester.addMappedFragment(0, 1, false, ELIGIBLE_BASE_QUALITY);
        tester.addMappedFragment(0, 1, true, DEFAULT_BASE_QUALITY); // duplicate!!!
        tester.addMappedFragment(1, 1, false, ELIGIBLE_BASE_QUALITY);
        tester.addMappedFragment(1, 1, true, DEFAULT_BASE_QUALITY); // duplicate!!!
        tester.runTest();
    }


    @Test
    public void testTwoGroupsOnDifferentChromosomesOfTwoMappedPairs() {
        final MarkDuplicatesSparkTester tester = getTester();
        tester.addMappedPair(0, 1, 100, false, false, ELIGIBLE_BASE_QUALITY);
        tester.addMappedPair(0, 1, 100, true, true, DEFAULT_BASE_QUALITY); // duplicate!!!
        tester.addMappedPair(1, 1, 100, false, false, ELIGIBLE_BASE_QUALITY);
        tester.addMappedPair(1, 1, 100, true, true, DEFAULT_BASE_QUALITY); // duplicate!!!
        tester.runTest();
    }

    @Test
    public void testTwoGroupsOnDifferentChromosomesOfThreeMappedPairs() {
        final MarkDuplicatesSparkTester tester = getTester();
        tester.addMappedPair(0, 1, 100, false, false, ELIGIBLE_BASE_QUALITY);
        tester.addMappedPair(0, 1, 100, true, true, DEFAULT_BASE_QUALITY); // duplicate!!!
        tester.addMappedPair(0, 1, 100, true, true, DEFAULT_BASE_QUALITY); // duplicate!!!
        tester.addMappedPair(1, 1, 100, false, false, ELIGIBLE_BASE_QUALITY);
        tester.addMappedPair(1, 1, 100, true, true, DEFAULT_BASE_QUALITY); // duplicate!!!
        tester.addMappedPair(1, 1, 100, true, true, DEFAULT_BASE_QUALITY); // duplicate!!!
        tester.runTest();
    }

    @Test
    public void testThreeGroupsOnDifferentChromosomesOfThreeMappedPairs() {
        final MarkDuplicatesSparkTester tester = getTester();
        tester.addMappedPair(0, 1, 100, false, false, ELIGIBLE_BASE_QUALITY);
        tester.addMappedPair(0, 1, 100, true, true, DEFAULT_BASE_QUALITY); // duplicate!!!
        tester.addMappedPair(0, 1, 100, true, true, DEFAULT_BASE_QUALITY); // duplicate!!!
        tester.addMappedPair(1, 1, 100, false, false, ELIGIBLE_BASE_QUALITY);
        tester.addMappedPair(1, 1, 100, true, true, DEFAULT_BASE_QUALITY); // duplicate!!!
        tester.addMappedPair(1, 1, 100, true, true, DEFAULT_BASE_QUALITY); // duplicate!!!
        tester.addMappedPair(2, 1, 100, false, false, ELIGIBLE_BASE_QUALITY);
        tester.addMappedPair(2, 1, 100, true, true, DEFAULT_BASE_QUALITY); // duplicate!!!
        tester.addMappedPair(2, 1, 100, true, true, DEFAULT_BASE_QUALITY); // duplicate!!!
        tester.runTest();
    }

    @Test
    public void testBulkFragmentsNoDuplicates() {
        final MarkDuplicatesSparkTester tester = getTester();
        tester.getSamRecordSetBuilder().setReadLength(100);
        for(int position = 1; position <= 10000; position += 1) {
            tester.addMappedFragment(0, position, false, "100M", DEFAULT_BASE_QUALITY);
        }
        tester.runTest();
    }

    @Test
    public void testBulkFragmentsWithDuplicates() {
        final MarkDuplicatesSparkTester tester = getTester();
        tester.getSamRecordSetBuilder().setReadLength(100);
        for(int position = 1; position <= 10000; position += 1) {
            tester.addMappedFragment(0, position, false, "100M", ELIGIBLE_BASE_QUALITY);
            tester.addMappedFragment(0, position, true, "100M", DEFAULT_BASE_QUALITY);
            tester.addMappedFragment(0, position, true, "100M", DEFAULT_BASE_QUALITY);
            tester.addMappedFragment(0, position, true, "100M", DEFAULT_BASE_QUALITY);
            tester.addMappedFragment(0, position, true, "100M", DEFAULT_BASE_QUALITY);
        }
        tester.runTest();
    }

    @Test
    public void testStackOverFlowPairSetSwap() {
        final MarkDuplicatesSparkTester tester = getTester();
        tester.getSamRecordSetBuilder().setReadLength(68);

        final File input = new File(TEST_DATA_DIR, "markDuplicatesWithMateCigar.pairSet.swap.sam");
        final SamReader reader = SamReaderFactory.makeDefault().open(input);
        tester.setHeader(reader.getFileHeader());
        for (final SAMRecord record : reader) {
            tester.addRecord(record);
        }
        CloserUtil.close(reader);
        tester.setExpectedOpticalDuplicate(1);
        tester.runTest();
    }

    @Test
    public void testSecondEndIsBeforeFirstInCoordinate() {
        final MarkDuplicatesSparkTester tester = getTester();
        tester.getSamRecordSetBuilder().setReadLength(68);
        tester.addMappedPair(0, 108855339, 108855323, false, false, "33S35M", "17S51M", true, true, false, DEFAULT_BASE_QUALITY); // +/-
        tester.runTest();
    }

    @Test
    public void testPathologicalOrderingAtTheSamePosition() {
        final MarkDuplicatesSparkTester tester = getTester();
        tester.getSamRecordSetBuilder().setReadLength(68);

        tester.setExpectedOpticalDuplicate(1);

        tester.addMatePair("RUNID:3:1:15013:113051", 0, 129384554, 129384554, false, false, false, false, "68M", "68M", false, false, false, false, false, ELIGIBLE_BASE_QUALITY);
        tester.addMatePair("RUNID:3:1:15029:113060", 0, 129384554, 129384554, false, false, true, true, "68M", "68M", false, false, false, false, false, DEFAULT_BASE_QUALITY);

        // Create the pathology
        try (final CloseableIterator<SAMRecord> iterator = tester.getRecordIterator()) {
            final int[] qualityOffset = {20, 30, 10, 40}; // creates an interesting pathological ordering
            int index = 0;
            while (iterator.hasNext()) {
                final SAMRecord record = iterator.next();
                final byte[] quals = new byte[record.getReadLength()];
                for (int i = 0; i < record.getReadLength(); i++) {
                    quals[i] = (byte) (qualityOffset[index] + 10);
                }
                record.setBaseQualities(quals);
                index++;
            }
        }

        // Run the test
        tester.runTest();
    }

    @Test
    public void testDifferentChromosomesInOppositeOrder() {
        final MarkDuplicatesSparkTester tester = getTester();
        tester.getSamRecordSetBuilder().setReadLength(101);
        tester.setExpectedOpticalDuplicate(1);
        tester.addMatePair("RUNID:6:101:17642:6835", 0, 1, 123989, 18281, false, false, true, true, "37S64M", "52M49S", false, false, false, false, false, DEFAULT_BASE_QUALITY, "1");
        tester.addMatePair("RUNID:6:101:17616:6888", 1, 0, 18281, 123989, false, false, false, false, "52M49S", "37S64M", false, false, false, false, false, ELIGIBLE_BASE_QUALITY, "1");
        tester.runTest();
    }

    @Test
    public void testOpticalDuplicateClustersAddingSecondEndFirstSameCoordinate() {
        final MarkDuplicatesSparkTester tester = getTester();
        tester.getSamRecordSetBuilder().setReadLength(68);
        tester.setExpectedOpticalDuplicate(1);
        tester.addMatePair("RUNID:1:1:15993:13361", 2, 41212324, 41212310, false, false, false, false, "33S35M", "19S49M", true, true, false, false, false, ELIGIBLE_BASE_QUALITY);
        tester.addMatePair("RUNID:1:1:16020:13352", 2, 41212324, 41212319, false, false, true, true, "33S35M", "28S40M", true, true, false, false, false, DEFAULT_BASE_QUALITY);
        tester.runTest();
    }

    @DataProvider(name = "secondarySupplementaryData")
    public Object[][] secondarySupplementaryData() {
        return new Object[][] {
                { true,  true , true},
                { true,  false, true},
                { false, true , true},
                { true,  true , false},
                { true,  false, false},
                { false, true , false}
        };
    }

    @Test(dataProvider = "secondarySupplementaryData")
    public void testTwoMappedPairsWithSupplementaryReads(final Boolean additionalFragSecondary, final Boolean additionalFragSupplementary, final Boolean fragLikeFirst) {
        final MarkDuplicatesSparkTester tester = getTester();
        tester.getSamRecordSetBuilder().setReadLength(68);
        tester.setExpectedOpticalDuplicate(1);
        tester.addMatePair("RUNID:1:1:15993:13361", 2, 41212324, 41212310, false, false, false, false, "33S35M", "19S49M", true, true, false, false, false, DEFAULT_BASE_QUALITY);
        tester.addMatePair("RUNID:1:1:16020:13352", 2, 41212324, 41212319, false, false, true, true, "33S35M", "28S40M", true, true, false, false, false, DEFAULT_BASE_QUALITY); // duplicate pair
        tester.addMappedFragment(fragLikeFirst ? "RUNID:1:1:15993:13361" : "RUNID:1:1:16020:13352", 1, 400, markSecondaryAndSupplementaryRecordsLikeTheCanonical() && !fragLikeFirst, null, null, additionalFragSecondary, additionalFragSupplementary, DEFAULT_BASE_QUALITY);
        tester.runTest();
    }

    @Test(dataProvider = "testDuplicateDetectionDataProviderWithMetrics")
    public void testDuplicateDetectionDataProviderWithMetrics(final File sam, final File expectedMetrics) throws IOException {
        final File outputDir = BaseTest.createTempDir("MarkDups");
        final File outputSam = new File(outputDir, "MarkDups" + ".sam");
        outputSam.deleteOnExit();
        final File metricsFile = new File(outputDir, "MarkDups" + ".duplicate_metrics");
        metricsFile.deleteOnExit();
        final CommandLineProgram markDuplicates = getCommandLineProgramInstance();
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("-I");
        args.add(sam.getAbsolutePath());
        args.add("-O");
        args.add(outputSam.getAbsolutePath());
        args.add("--"+StandardArgumentDefinitions.METRICS_FILE_LONG_NAME);
        args.add(metricsFile.getAbsolutePath());
        markDuplicates.instanceMain(args.getArgsArray());
        IntegrationTestSpec.assertEqualTextFiles(metricsFile, expectedMetrics, "#"); //this compares the values but not headers

        //Note: headers need to be compares not by exact values because they include times and class names
        final List<String> lines = FileUtils.readLines(metricsFile, StandardCharsets.UTF_8);
        Assert.assertTrue(lines.get(0).startsWith("##"), lines.get(0));
        Assert.assertTrue(lines.get(1).startsWith("#"), lines.get(1));
        Assert.assertTrue(lines.get(1).toLowerCase().contains("--input"), lines.get(1));  //Note: lowercase because picard uses INPUT and GATK uses input for full name
        Assert.assertTrue(lines.get(2).startsWith("##"), lines.get(2));
        Assert.assertTrue(lines.get(3).startsWith("# Started on:"), lines.get(3));
        Assert.assertTrue(lines.get(4).trim().isEmpty());
        Assert.assertTrue(lines.get(4).trim().isEmpty());
        Assert.assertTrue(lines.get(5).startsWith("## METRICS CLASS"), lines.get(5));
    }

    @DataProvider(name="testDuplicateDetectionDataProviderWithMetrics")
    public Object[][] testDuplicateDetectionDataProviderWithMetrics() {
        //Note: the expected metrics files were created using picard 1.130
        return new Object[][] {
                {new File(TEST_DATA_DIR, "example.chr1.1-1K.unmarkedDups.bam"), new File(TEST_DATA_DIR,"expected.picard-2.15.0.sorted.chr1.1-1K.unmarkedDups.markDuplicate.metrics")},
                {new File(TEST_DATA_DIR, "optical_dupes.bam"), new File(TEST_DATA_DIR,"expected.picard-2.15.0.sorted.optical_dupes.markDuplicate.metrics")},
                {new File(TEST_DATA_DIR, "inputSingleLibrarySolexa16404.bam"), new File(TEST_DATA_DIR,"expected.picard-2.15.0.sorted.inputSingleLibrarySolexa16404.metrics")},
        };
    }

    @DataProvider(name="testMDdata")
    public Object[][] testMDdata() {
        //Note: the expected metrics files were created using picard 1.130

        //On this file, the spark version used to create the wrong order of reads.
        return new Object[][] {
                {new File(TEST_DATA_DIR, "mdOrderBug.bam"), new File(TEST_DATA_DIR,"expected.mdOrderBug.bam")},  //fixed in https://github.com/broadinstitute/gatk/pull/1197
                {new File(TEST_DATA_DIR, "mdOrderBug2.bam"), new File(TEST_DATA_DIR,"expected.mdOrderBug2.bam")},
        };
    }

    @Test(dataProvider = "testMDdata")
    public void testMDOrder(final File input, final File expectedOutput) throws Exception {
        // This method is overridden in MarkDuplicatesSparkIntegrationTest to provide a --num-reducers argument
        testMDOrderImpl(input, expectedOutput, "");
    }

    @Test
    public void testTwoMappedPairsWithSoftClippingFirstOfPairOnly() {
        final MarkDuplicatesSparkTester tester = getTester();
        // NB: no duplicates
        // 5'1: 2, 5'2:46+73M=118
        // 5'1: 2, 5'2:51+68M=118
        tester.getSamRecordSetBuilder().setReadLength(76);
        tester.addMappedPair(0, 12, 46, false, false, "6S42M28S", "3S73M", true, 50); // only add the first one
        // NB: this next record should not be a duplicate in MarkDuplicatesGATK
        tester.addMappedPair(0, 12, 51, false, false, "6S42M28S", "8S68M", true, 50); // only add the first one
        tester.runTest();
    }

    protected void testMDOrderImpl(final File input, final File expectedOutput, final String extraArgs) throws Exception {
        final File metricsFile = createTempFile("markdups_metrics", ".txt");
        final File outputFile = createTempFile("markdups", ".bam");

        ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("-"+ StandardArgumentDefinitions.INPUT_SHORT_NAME);
        args.add(input.getPath());
        args.add("-"+StandardArgumentDefinitions.OUTPUT_SHORT_NAME);
        args.add(outputFile.getAbsolutePath());
        args.add("--"+StandardArgumentDefinitions.METRICS_FILE_LONG_NAME);
        args.add(metricsFile.getAbsolutePath());
        args.add(extraArgs);

        final CommandLineProgram markDuplicates = getCommandLineProgramInstance();

        markDuplicates.instanceMain(args.getArgsArray());
        SamAssertionUtils.assertEqualBamFiles(outputFile, expectedOutput, false, ValidationStringency.SILENT);
    }


    @Test
    public void testNonExistantReadGroupInRead() {
        final MarkDuplicatesSparkTester tester = new MarkDuplicatesSparkTester(true);
        tester.addMatePair("RUNID:7:1203:2886:82292",  19, 19, 485253, 485253, false, false, true, true, "42M59S", "59S42M", true, false, false, false, false, DEFAULT_BASE_QUALITY, "NotADuplicateGroup");
        try {
            tester.runTest();
            Assert.fail("Should have thrown an exception");
        } catch (Exception e){
           Assert.assertTrue(e instanceof SparkException);
           Assert.assertTrue(e.getCause() instanceof UserException.HeaderMissingReadGroup);
        }
    }

    @Test
    public void testNoReadGroupInRead() {
        final MarkDuplicatesSparkTester tester = new MarkDuplicatesSparkTester(true);
        tester.addMatePair("RUNID:7:1203:2886:82292",  19, 19, 485253, 485253, false, false, true, true, "42M59S", "59S42M", true, false, false, false, false, DEFAULT_BASE_QUALITY, null);

        try {
            tester.runTest();
            Assert.fail("Should have thrown an exception");
        } catch (Exception e){
            Assert.assertTrue(e instanceof SparkException);
            Assert.assertTrue(e.getCause() instanceof UserException.ReadMissingReadGroup);
        }
    }

    @Test(dataProvider = "readNameData")
    public void testOpticalDuplicatesTiebrokenByPhysicalLocationNotStartPosition(final String readName1, final String readName2) {
        // This tests the readname based tiebreaking code in mark duplicates. Since it's ambiguous which read should be marked
        // as duplicate or not if scores match we break ties by evaluating the readname for consistencies sake.

        final ReadNameParser parser = new ReadNameParser();

        final PhysicalLocationInt position1 = new PhysicalLocationInt();
        final PhysicalLocationInt position2 = new PhysicalLocationInt();

        parser.addLocationInformation(readName1, position1);
        parser.addLocationInformation(readName2, position2);

        final MarkDuplicatesSparkTester tester = getTester();
        tester.getSamRecordSetBuilder().setReadLength(101);
        tester.setExpectedOpticalDuplicate(0);

        int compare = position1.tile - position2.tile;
        if (compare == 0) {
            compare = (short)position1.x - (short)position2.x;
        }

        if (compare == 0) {
            compare = (short)position1.y - (short)position2.y;
        }

        final boolean isDuplicate = compare < 0;

        // NOTE these reads are offset slightly but should have the same unclipped start postitions
        tester.addMatePair(readName1, 1,2, 46,  false, false, !isDuplicate, !isDuplicate, "6S42M28S", "3S68M",  true, false, false, false, false, DEFAULT_BASE_QUALITY);
        tester.addMatePair(readName2, 1,2, 51, false, false, isDuplicate, isDuplicate, "6S42M28S", "8S68M", true, false, false, false, false, DEFAULT_BASE_QUALITY);
        tester.runTest();
    }

    @DataProvider
    public Object[][] readNameData(){
        return new Object[][]{
                {"RUNID:7:1203:2886:82292", "RUNID:7:1205:3886:16834"},

                {"RUNID:7:1203:2886:16756", "RUNID:7:1205:3886:16756"},
                {"RUNID:7:1204:2886:16756", "RUNID:7:1205:3886:16756"},
                {"RUNID:7:1205:2886:16756", "RUNID:7:1205:3886:16756"},
                {"RUNID:7:1206:2886:16756", "RUNID:7:1205:3886:16756"},
                {"RUNID:7:1207:2886:16756", "RUNID:7:1205:3886:16756"},

                {"RUNID:7:1203:2886:16756", "RUNID:7:1203:4886:26756"},
                {"RUNID:7:1203:3886:16756", "RUNID:7:1203:4886:26756"},
                {"RUNID:7:1203:4886:16756", "RUNID:7:1203:4886:26756"},
                {"RUNID:7:1203:5886:16756", "RUNID:7:1203:4886:26756"},
                {"RUNID:7:1203:6886:16756", "RUNID:7:1203:4886:26756"},

                {"RUNID:7:1203:2886:34756", "RUNID:7:1203:2886:36756"},
                {"RUNID:7:1203:2886:35756", "RUNID:7:1203:2886:36756"},
                {"RUNID:7:1203:2886:37756", "RUNID:7:1203:2886:36756"},
                {"RUNID:7:1203:2886:38756", "RUNID:7:1203:2886:36756"},

                //Added a test for tiebreaking accounting for the short casting done in picard
                {"HK3T5CCXX160204:3:1112:11586:37067", "HK3T5CCXX160204:3:1112:11586:32144"}
        };
    }

    @Test(dataProvider = "readNameData")
    public void testOpticalDuplicateClusterSamePositionNoOpticalDuplicates(final String readName1, final String readName2) {
        // This tests the readname based tiebreaking code in mark duplicates. Since it's ambiguous which read should be marked
        // as duplicate or not if scores match we break ties by evaluating the readname for consistencies sake.

        final ReadNameParser parser = new ReadNameParser();

        final PhysicalLocationInt position1 = new PhysicalLocationInt();
        final PhysicalLocationInt position2 = new PhysicalLocationInt();

        parser.addLocationInformation(readName1, position1);
        parser.addLocationInformation(readName2, position2);

        final MarkDuplicatesSparkTester tester = getTester();
        tester.getSamRecordSetBuilder().setReadLength(101);
        tester.setExpectedOpticalDuplicate(0);

        int compare = position1.tile - position2.tile;
        if (compare == 0) {
            compare = (short)position1.x - (short)position2.x;
        }

        if (compare == 0) {
            compare = (short)position1.y - (short)position2.y;
        }

        final boolean isDuplicate = compare < 0;

        tester.addMatePair(readName1, 1,485253, 485253, false, false, !isDuplicate, !isDuplicate, "42M59S", "59S42M", false, true, false, false, false, DEFAULT_BASE_QUALITY);
        tester.addMatePair(readName2, 1,485253, 485253, false, false, isDuplicate, isDuplicate, "59S42M", "42M59S", true, false, false, false, false, DEFAULT_BASE_QUALITY);

        tester.runTest();
    }
}
