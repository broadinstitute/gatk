package org.broadinstitute.hellbender.tools.picard.sam;

import htsjdk.samtools.*;
import htsjdk.samtools.util.CloserUtil;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.test.SamAssertionUtils;
import org.broadinstitute.hellbender.utils.read.mergealignment.BestMapqPrimaryAlignmentSelectionStrategy;
import org.broadinstitute.hellbender.utils.read.mergealignment.SamAlignmentMerger;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.*;

/**
 *  Test for the MergeBamAlignment class
 *
 * @author ktibbett@broadinstitute.org
 */
public final class MergeBamAlignmentIntegrationTest extends CommandLineProgramTest {

    private static final File TEST_DATA_DIR = new File(getTestDataDir(), "picard/sam/MergeBamAlignment");

    private static final File unmappedCram = new File(TEST_DATA_DIR, "unmapped.cram");
    private static final File alignedCram = new File(TEST_DATA_DIR, "aligned.cram");
    private static final File unmappedBam = new File(TEST_DATA_DIR, "unmapped.sam");
    private static final File alignedBam = new File(TEST_DATA_DIR, "aligned.sam");
    private static final File oneHalfAlignedBam = new File(TEST_DATA_DIR, "onehalfaligned.sam");
    private static final File otherHalfAlignedBam = new File(TEST_DATA_DIR, "otherhalfaligned.sam");
    private static final File mergingUnmappedBam = new File(TEST_DATA_DIR, "merging.unmapped.sam");
    private static final File firstReadAlignedBam = new File(TEST_DATA_DIR, "allread1.trimmed.aligned.sam");
    private static final File secondReadAlignedBam = new File(TEST_DATA_DIR, "allread2.trimmed.aligned.sam");
    private static final File firstReadAlignedBam_firstHalf = new File(TEST_DATA_DIR, "firsthalf.read1.trimmed.aligned.sam");
    private static final File firstReadAlignedBam_secondHalf = new File(TEST_DATA_DIR, "secondhalf.read1.trimmed.aligned.sam");
    private static final File secondReadAlignedBam_firstHalf = new File(TEST_DATA_DIR, "firsthalf.read2.trimmed.aligned.sam");
    private static final File secondReadAlignedBam_secondHalf = new File(TEST_DATA_DIR, "secondhalf.read2.trimmed.aligned.sam");
    private static final File supplementalReadAlignedBam = new File(TEST_DATA_DIR, "aligned.supplement.sam");
    private static final File badorderUnmappedBam = new File(TEST_DATA_DIR, "unmapped.badorder.sam");
    private static final File badorderAlignedBam = new File(TEST_DATA_DIR, "aligned.badorder.sam");
    private static final File alignedQuerynameSortedBam = new File(TEST_DATA_DIR, "aligned_queryname_sorted.sam");
    private static final File multiHitUnmappedBam = new File(TEST_DATA_DIR, "multihit.unmapped.sam");
    private static final File multiHitAlignedBam = new File(TEST_DATA_DIR, "multihit.aligned.sam");
    private static final File multiHitFilterUnmappedBam = new File(TEST_DATA_DIR, "multihit.filter.unmapped.sam");
    private static final File multiHitFilterFragmentUnmappedBam = new File(TEST_DATA_DIR, "multihit.filter.fragment.unmapped.sam");
    private static final File clipTestUnmappedBam = new File(TEST_DATA_DIR, "cliptest.unmapped.sam");
    private static final File clipTestAlignedBam = new File(TEST_DATA_DIR, "cliptest.aligned.sam");
    private static final File clippedReference = new File(TEST_DATA_DIR, "cliptest.fasta");

    private static final File fasta = new File(TEST_DATA_DIR, "merger.fasta");
    private static final String bigSequenceName = "chr7"; // The longest sequence in merger.fasta
    private static final File sequenceDict = new File(TEST_DATA_DIR, "merger.dict");

    // For EarliestFragment tests, tag placed on the alignments which are expected to be marked as primary.
    private static final String ONE_OF_THE_BEST_TAG = "YB";

    public String getTestedClassName() {
        return MergeBamAlignment.class.getSimpleName();
    }

    @Test
    public void testMergerWithSupplemental() throws Exception {
        final File outputWithSupplemental = BaseTest.createTempFile("mergeWithSupplementalTest", ".sam");
        System.out.println(outputWithSupplemental.getAbsolutePath());
        // outputWithSupplemental.deleteOnExit();

        doMergeAlignment(unmappedBam,
                Collections.singletonList(supplementalReadAlignedBam),
                null, null, null, null,
                false, true, false, 1,
                "0", "1.0", "align!", "myAligner",
                fasta, outputWithSupplemental,
                SamPairUtil.PairOrientation.FR, null, null, null
        );

        final SamReader result = SamReaderFactory.makeDefault().open(outputWithSupplemental);

        final List<Integer> clipAdapterFlags = new ArrayList<>(Arrays.asList(99, 2147, 147, 2195));
        final List<Integer> foundClipAdapterFlags = new ArrayList<>();

        for (final SAMRecord sam : result) {
            if (sam.getReadName().equals("both_reads_align_clip_adapter")) {
                foundClipAdapterFlags.add(sam.getFlags());
            }

            // This tests that we clip both (a) when the adapter is marked in the unmapped BAM file and
            // (b) when the insert size is less than the read length
            if (sam.getReadName().equals("both_reads_align_clip_marked")) {
                Assert.assertEquals(sam.getReferenceName(), "chr7");
                if (sam.getReadNegativeStrandFlag()) {
                    Assert.assertEquals(sam.getCigarString(), "5S96M", "Incorrect CIGAR string for " + sam.getReadName());
                } else {
                    Assert.assertEquals(sam.getCigarString(), "96M5S", "Incorrect CIGAR string for " + sam.getReadName());
                }
            }
            else if (sam.getReadName().equals("both_reads_align_clip_adapter")) {
                Assert.assertEquals(sam.getReferenceName(), "chr7");
                if (!sam.getSupplementaryAlignmentFlag()) {
                    if (sam.getReadNegativeStrandFlag()) Assert.assertEquals(sam.getCigarString(), "5S96M", "Incorrect CIGAR string for " + sam.getReadName());
                    else Assert.assertEquals(sam.getCigarString(), "96M5S", "Incorrect CIGAR string for " + sam.getReadName());
                }
            }
            // This tests that we DON'T clip when we run off the end if there are equal to or more than
            // MIN_ADAPTER_BASES hanging off the end
            else if (sam.getReadName().equals("both_reads_align_min_adapter_bases_exceeded")) {
                Assert.assertEquals(sam.getReferenceName(), "chr7");
                Assert.assertTrue(!sam.getCigarString().contains("S"),
                        "Read was clipped when it should not be.");
            } else if (sam.getReadName().equals("neither_read_aligns_or_present")) {
                Assert.assertTrue(sam.getReadUnmappedFlag(), "Read should be unmapped but isn't");
            }
            // Two pairs in which only the first read should align
            else if (sam.getReadName().equals("both_reads_present_only_first_aligns") ||
                    sam.getReadName().equals("read_2_too_many_gaps")) {
                if (sam.getFirstOfPairFlag()) {
                    Assert.assertEquals(sam.getReferenceName(), "chr7", "Read should be mapped but isn't");
                } else {
                    Assert.assertTrue(sam.getReadUnmappedFlag(), "Read should not be mapped but is");
                }
            } else {
                throw new Exception("Unexpected read name: " + sam.getReadName());
            }
        }

        // Make sure that we have the appropriate primary and supplementary reads in the new file
        Assert.assertEquals(foundClipAdapterFlags.size(), clipAdapterFlags.size());
        Collections.sort(clipAdapterFlags);
        Collections.sort(foundClipAdapterFlags);
        for (int i = 0; i < clipAdapterFlags.size(); i++) {
            Assert.assertEquals(clipAdapterFlags.get(i), foundClipAdapterFlags.get(i));
        }

    }

    private void doMergerTest(final File unmappedInputFile,final File alignedInputFile, final String outputExtension) throws Exception {
        final File output = BaseTest.createTempFile("mergeTest", outputExtension);

        doMergeAlignment(unmappedInputFile, Collections.singletonList(alignedInputFile),
                null, null, null, null,
                false, true, false, 1,
                "0", "1.0", "align!", "myAligner",
                fasta, output,
                SamPairUtil.PairOrientation.FR, null, null, null
        );

        SamAssertionUtils.assertCRAMContentsIfCRAM(output);

        SamReader result = SamReaderFactory.makeDefault().referenceSequence(fasta).open(output);
        Assert.assertEquals(result.getFileHeader().getSequenceDictionary().getSequences().size(), 8,
                "Number of sequences did not match");
        SAMProgramRecord pg = result.getFileHeader().getProgramRecords().get(0);
        Assert.assertEquals(pg.getProgramGroupId(), "0");
        Assert.assertEquals(pg.getProgramVersion(), "1.0");
        Assert.assertEquals(pg.getCommandLine(), "align!");
        Assert.assertEquals(pg.getProgramName(), "myAligner");

        final SAMReadGroupRecord rg = result.getFileHeader().getReadGroups().get(0);
        Assert.assertEquals(rg.getReadGroupId(), "0");
        Assert.assertEquals(rg.getSample(), "Hi,Mom!");

        Assert.assertEquals(result.getFileHeader().getSortOrder(), SAMFileHeader.SortOrder.coordinate);

        for (final SAMRecord sam : result) {
            // This tests that we clip both (a) when the adapter is marked in the unmapped BAM file and
            // (b) when the insert size is less than the read length
            if (sam.getReadName().equals("both_reads_align_clip_adapter") ||
                    sam.getReadName().equals("both_reads_align_clip_marked")) {
                Assert.assertEquals(sam.getReferenceName(), "chr7");
                if (sam.getReadNegativeStrandFlag()) {
                    Assert.assertEquals(sam.getCigarString(), "5S96M", "Incorrect CIGAR string for " +
                            sam.getReadName());
                } else {
                    Assert.assertEquals(sam.getCigarString(), "96M5S", "Incorrect CIGAR string for " +
                            sam.getReadName());
                }
            }
            // This tests that we DON'T clip when we run off the end if there are equal to or more than
            // MIN_ADAPTER_BASES hanging off the end
            else if (sam.getReadName().equals("both_reads_align_min_adapter_bases_exceeded")) {
                Assert.assertEquals(sam.getReferenceName(), "chr7");
                Assert.assertTrue(!sam.getCigarString().contains("S"),
                        "Read was clipped when it should not be.");
            } else if (sam.getReadName().equals("neither_read_aligns_or_present")) {
                Assert.assertTrue(sam.getReadUnmappedFlag(), "Read should be unmapped but isn't");
            }
            // Two pairs in which only the first read should align
            else if (sam.getReadName().equals("both_reads_present_only_first_aligns") ||
                    sam.getReadName().equals("read_2_too_many_gaps")) {
                if (sam.getFirstOfPairFlag()) {
                    Assert.assertEquals(sam.getReferenceName(), "chr7", "Read should be mapped but isn't");
                } else {
                    Assert.assertTrue(sam.getReadUnmappedFlag(), "Read should not be mapped but is");
                }
            } else {
                throw new Exception("Unexpected read name: " + sam.getReadName());
            }

        }
    }

    @Test
    public void testCRAMMerger() throws Exception {
        doMergerTest(unmappedCram, alignedCram, ".cram");
    }

    @Test
    public void testBAMMerger() throws Exception {
        doMergerTest(unmappedBam, alignedBam, ".sam");
    }

    @Test
    public void testCRAMAndBAMMerger() throws Exception {
        doMergerTest(unmappedCram, alignedBam, ".sam");
    }

    @Test
    public void testMergerFromMultipleFiles() throws Exception {
        final File output = BaseTest.createTempFile("mergeTest", ".sam");

        doMergeAlignment(unmappedBam, Arrays.asList(oneHalfAlignedBam, otherHalfAlignedBam),
                null, null, null, null,
                false, true, false, 1,
                "0", "1.0", "align!", "myAligner",
                fasta, output,
                SamPairUtil.PairOrientation.FR, null, null, null
        );

        final SamReader result = SamReaderFactory.makeDefault().open(output);

        for (final SAMRecord sam : result) {
            // This tests that we clip both (a) when the adapter is marked in the unmapped BAM file and
            // (b) when the insert size is less than the read length
            if (sam.getReadName().equals("both_reads_align_clip_adapter") ||
                    sam.getReadName().equals("both_reads_align_clip_marked")) {
                Assert.assertEquals(sam.getReferenceName(), "chr7");
                if (sam.getReadNegativeStrandFlag()) {
                    Assert.assertEquals(sam.getCigarString(), "5S96M", "Incorrect CIGAR string for " +
                            sam.getReadName());
                } else {
                    Assert.assertEquals(sam.getCigarString(), "96M5S", "Incorrect CIGAR string for " +
                            sam.getReadName());
                }
            }
            // This tests that we DON'T clip when we run off the end if there are equal to or more than
            // MIN_ADAPTER_BASES hanging off the end
            else if (sam.getReadName().equals("both_reads_align_min_adapter_bases_exceeded")) {
                Assert.assertEquals(sam.getReferenceName(), "chr7");
                Assert.assertTrue(!sam.getCigarString().contains("S"),
                        "Read was clipped when it should not be.");
            } else if (sam.getReadName().equals("neither_read_aligns_or_present")) {
                Assert.assertTrue(sam.getReadUnmappedFlag(), "Read should be unmapped but isn't");
            }
            // Two pairs in which only the first read should align
            else if (sam.getReadName().equals("both_reads_present_only_first_aligns") ||
                    sam.getReadName().equals("read_2_too_many_gaps")) {
                if (sam.getFirstOfPairFlag()) {
                    Assert.assertEquals(sam.getReferenceName(), "chr7", "Read should be mapped but isn't");
                } else {
                    Assert.assertTrue(sam.getReadUnmappedFlag(), "Read should not be mapped but is");
                }
            } else {
                throw new Exception("Unexpected read name: " + sam.getReadName());
            }

        }

    }

    @Test(dataProvider="data")
    public void testSortingOnSamAlignmentMerger(final File unmapped, final File aligned, final boolean sorted, final String testName)
            throws IOException {

        final File target = BaseTest.createTempFile("target", "bam");
        final SamAlignmentMerger merger = new SamAlignmentMerger(unmapped,  target, fasta, null, true, false,
                false, Arrays.asList(aligned), 1, null, null, null, null, null, null,
                Arrays.asList(SamPairUtil.PairOrientation.FR), SAMFileHeader.SortOrder.coordinate,
                new BestMapqPrimaryAlignmentSelectionStrategy(), false);

        merger.mergeAlignment(Defaults.REFERENCE_FASTA);
        Assert.assertEquals(sorted, !merger.getForceSort());
        final SAMRecordIterator it = SamReaderFactory.makeDefault().open(target).iterator();
        int aln = 0;
        while (it.hasNext()) {
            final SAMRecord rec = it.next();
            if (!rec.getReadUnmappedFlag()) {
                aln++;
            }
        }
        Assert.assertEquals(aln, 6, "Incorrect number of aligned reads in merged BAM file");
    }

    @DataProvider(name="data")
    public Object[][] getDataForSortingTest() {
        return new Object[][] {
                {unmappedBam, alignedQuerynameSortedBam, true, "Basic test with pre-sorted alignment"},
                {unmappedBam, alignedBam, false, "Basic test with unsorted alignment"}

        };
    }

    /**
     * Minimal test of merging data from separate read 1 and read 2 alignments
     */
    @Test(dataProvider="separateTrimmed")
    public void testMergingFromSeparatedReadTrimmedAlignments(final File unmapped, final List<File> r1Align, final List<File> r2Align, final int r1Trim, final int r2Trim, final String testName) throws Exception {
        final File output = BaseTest.createTempFile("mergeMultipleAlignmentsTest", ".sam");

        doMergeAlignment(unmapped, null, r1Align, r2Align, r1Trim, r2Trim,
                false, true, false, 1,
                "0", "1.0", "align!", "myAligner",
                fasta, output,
                SamPairUtil.PairOrientation.FR, null, null, null
        );

        SamReaderFactory factory = SamReaderFactory.makeDefault();
        final SamReader result = factory.open(output);
        for (final SAMRecord sam : result) {
            // Get the alignment record
            final List<File> rFiles = sam.getFirstOfPairFlag() ? r1Align : r2Align;
            SAMRecord alignment = null;
            for (final File f : rFiles) {
                for (final SAMRecord tmp : factory.open(f)) {
                    if (tmp.getReadName().equals(sam.getReadName())) {
                        alignment = tmp;
                        break;
                    }
                }
                if (alignment != null) break;
            }

            final int trim = sam.getFirstOfPairFlag() ? r1Trim : r2Trim;
            final int notWrittenToFastq = sam.getReadLength() - (trim + alignment.getReadLength());
            final int beginning = sam.getReadNegativeStrandFlag() ? notWrittenToFastq : trim;
            final int end = sam.getReadNegativeStrandFlag() ? trim : notWrittenToFastq;

            if (!sam.getReadUnmappedFlag()) {
                final CigarElement firstMergedCigarElement = sam.getCigar().getCigarElement(0);
                final CigarElement lastMergedCigarElement = sam.getCigar().getCigarElement(sam.getCigar().numCigarElements() - 1);
                final CigarElement firstAlignedCigarElement = alignment.getCigar().getCigarElement(0);
                final CigarElement lastAlignedCigarElement = alignment.getCigar().getCigarElement(alignment.getCigar().numCigarElements() - 1);

                if (beginning > 0) {
                    Assert.assertEquals(firstMergedCigarElement.getOperator(), CigarOperator.S, "First element is not a soft clip");
                    Assert.assertEquals(firstMergedCigarElement.getLength(), beginning + ((firstAlignedCigarElement.getOperator() == CigarOperator.S) ? firstAlignedCigarElement.getLength() : 0));
                }
                if (end > 0) {
                    Assert.assertEquals(lastMergedCigarElement.getOperator(), CigarOperator.S, "Last element is not a soft clip");
                    Assert.assertEquals(lastMergedCigarElement.getLength(), end + ((lastAlignedCigarElement.getOperator() == CigarOperator.S) ? lastAlignedCigarElement.getLength() : 0));
                }
            }
        }

    }

    @DataProvider(name="separateTrimmed")
    public Object[][] getDataForMergingTest() {
        return new Object[][] {
                // Tests using one file for each read, queryname sorted, trimmed, not all bases written
                {mergingUnmappedBam, Arrays.asList(firstReadAlignedBam), Arrays.asList(secondReadAlignedBam), 17, 20, "one file per read"},
                // Tests using multiple files for each read, coordinate sorted, trimmed, not all bases written
                {mergingUnmappedBam, Arrays.asList(firstReadAlignedBam_firstHalf, firstReadAlignedBam_secondHalf),
                        Arrays.asList(secondReadAlignedBam_firstHalf, secondReadAlignedBam_secondHalf), 17, 20, "two files per read"}
        };
    }

    @Test(expectedExceptions = {IllegalStateException.class, GATKException.class})
    public void testOldQuerynameSortFails() throws IOException {

        final File merged = BaseTest.createTempFile("merged", BamFileIoUtils.BAM_FILE_EXTENSION);

        doMergeAlignment(badorderUnmappedBam, Collections.singletonList(badorderAlignedBam),
                null, null, null, null,
                false, true, false, 1,
                "0", "1.0", "align!", "myAligner",
                fasta, merged,
                SamPairUtil.PairOrientation.FR, null,
                null, null
        );
        Assert.fail("Merger should have failed because unmapped reads are not in queryname order but didn't");
    }

    @Test
    public void testMultiHit() throws IOException {
        final File merged = BaseTest.createTempFile("merged", ".sam");

        doMergeAlignment(multiHitUnmappedBam, Collections.singletonList(multiHitAlignedBam),
                null, null, null, null,
                false, true, false, 1,
                "0", "1.0", "align!", "myAligner",
                fasta, merged,
                null, null, null, null
        );

        // Iterate over the merged output and gather some statistics
        final Map<String, AlignmentAccumulator> accumulatorMap = new HashMap<>();

        final SamReader reader = SamReaderFactory.makeDefault().open(merged);
        for (final SAMRecord rec : reader) {
            final String readName;
            if (!rec.getReadPairedFlag()) readName = rec.getReadName();
            else if (rec.getFirstOfPairFlag()) readName = rec.getReadName() + "/1";
            else readName = rec.getReadName() + "/2";
            AlignmentAccumulator acc = accumulatorMap.get(readName);
            if (acc == null) {
                acc = new AlignmentAccumulator();
                accumulatorMap.put(readName, acc);
            }
            if (rec.getReadUnmappedFlag()) {
                Assert.assertFalse(acc.seenUnaligned, readName);
                acc.seenUnaligned = true;
            } else {
                ++acc.numAlignments;
                if (!rec.getNotPrimaryAlignmentFlag()) {
                    Assert.assertNull(acc.primaryAlignmentSequence, readName);
                    acc.primaryAlignmentSequence = rec.getReferenceName();
                }
            }
        }

        // Set up expected results
        final Map<String, AlignmentAccumulator> expectedMap = new HashMap<>();
        expectedMap.put("frag_hit", new AlignmentAccumulator(false, 1, "chr1"));
        expectedMap.put("frag_multihit", new AlignmentAccumulator(false, 3, "chr2"));
        expectedMap.put("frag_no_hit", new AlignmentAccumulator(true, 0, null));
        expectedMap.put("pair_both_hit/1", new AlignmentAccumulator(false, 1, "chr7"));
        expectedMap.put("pair_both_hit/2", new AlignmentAccumulator(false, 1, "chr7"));
        expectedMap.put("pair_both_multihit/1", new AlignmentAccumulator(false, 3, "chr8"));
        expectedMap.put("pair_both_multihit/2", new AlignmentAccumulator(false, 3, "chr8"));
        expectedMap.put("pair_first_hit/1", new AlignmentAccumulator(false, 1, "chr1"));
        expectedMap.put("pair_first_hit/2", new AlignmentAccumulator(true, 0, null));
        expectedMap.put("pair_first_multihit/1", new AlignmentAccumulator(false, 3, "chr1"));
        expectedMap.put("pair_first_multihit/2", new AlignmentAccumulator(true, 0, null));
        expectedMap.put("pair_no_hit/1", new AlignmentAccumulator(true, 0, null));
        expectedMap.put("pair_no_hit/2", new AlignmentAccumulator(true, 0, null));
        expectedMap.put("pair_second_hit/1", new AlignmentAccumulator(true, 0, null));
        expectedMap.put("pair_second_hit/2", new AlignmentAccumulator(false, 1, "chr1"));
        expectedMap.put("pair_second_multihit/1", new AlignmentAccumulator(true, 0, null));
        expectedMap.put("pair_second_multihit/2", new AlignmentAccumulator(false, 3, "chr3"));

        Assert.assertEquals(accumulatorMap.size(), expectedMap.size());

        for (final Map.Entry<String, AlignmentAccumulator> entry : expectedMap.entrySet()) {
            final AlignmentAccumulator actual = accumulatorMap.get(entry.getKey());
            Assert.assertEquals(actual, entry.getValue(), entry.getKey());
        }

        SamAssertionUtils.assertSamValid(merged);
    }

    private static class AlignmentAccumulator {
        boolean seenUnaligned = false;
        int numAlignments = 0;
        String primaryAlignmentSequence = null;

        private AlignmentAccumulator() {}

        private AlignmentAccumulator(final boolean seenUnaligned, final int numAlignments, final String primaryAlignmentSequence) {
            this.seenUnaligned = seenUnaligned;
            this.numAlignments = numAlignments;
            this.primaryAlignmentSequence = primaryAlignmentSequence;
        }

        public boolean equals(final Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            final AlignmentAccumulator that = (AlignmentAccumulator) o;

            if (numAlignments != that.numAlignments) return false;
            if (seenUnaligned != that.seenUnaligned) return false;
            if (primaryAlignmentSequence != null ? !primaryAlignmentSequence.equals(that.primaryAlignmentSequence) : that.primaryAlignmentSequence != null)
                return false;

            return true;
        }

        @Override
        public int hashCode() {
            int result = (seenUnaligned ? 1 : 0);
            result = 31 * result + numAlignments;
            result = 31 * result + (primaryAlignmentSequence != null ? primaryAlignmentSequence.hashCode() : 0);
            return result;
        }

        @Override
        public String toString() {
            return "AlignmentAccumulator{" +
                    "seenUnaligned=" + seenUnaligned +
                    ", numAlignments=" + numAlignments +
                    ", primaryAlignmentSequence='" + primaryAlignmentSequence + '\'' +
                    '}';
        }
    }

    /**
     * Read a single paired-end read from a file, and create one or more aligned records for the read pair based on
     * the lists, merge with the original paired-end read, and assert expected results.
     * @param description
     * @param firstOfPair List that describes the aligned SAMRecords to create for the first end.
     * @param secondOfPair List that describes the aligned SAMRecords to create for the second end.
     * @param expectedPrimaryHitIndex Expected value for the HI tag in the primary alignment in the merged output.
     * @param expectedNumFirst Expected number of first-of-pair SAMRecords in the merged output.
     * @param expectedNumSecond Expected number of second-of-pair SAMRecords in the merged output.
     * @param expectedPrimaryMapq Sum of mapqs of both ends of primary alignment in the merged output.
     * @throws Exception
     */
    @Test(dataProvider = "testPairedMultiHitWithFilteringTestCases")
    public void testPairedMultiHitWithFiltering(final String description, final List<HitSpec> firstOfPair, final List<HitSpec> secondOfPair,
                                                final Integer expectedPrimaryHitIndex, final int expectedNumFirst,
                                                final int expectedNumSecond, final int expectedPrimaryMapq) throws Exception {

        SAMRecord firstUnmappedRec = null;
        SAMRecord secondUnmappedRec = null;

        // Create the aligned file by copying bases, quals, readname from the unmapped read, and conforming to each HitSpec.
        try (final SAMRecordIterator unmappedSamFileIterator = SamReaderFactory.makeDefault().open(multiHitFilterUnmappedBam).iterator()) {
            firstUnmappedRec = unmappedSamFileIterator.next();
            secondUnmappedRec = unmappedSamFileIterator.next();
        }
        final File alignedSam = BaseTest.createTempFile("aligned.", ".sam");
        final SAMFileHeader alignedHeader = new SAMFileHeader();
        alignedHeader.setSortOrder(SAMFileHeader.SortOrder.queryname);
        alignedHeader.setSequenceDictionary(SamReaderFactory.makeDefault().getFileHeader(sequenceDict).getSequenceDictionary());
        try (final SAMFileWriter alignedWriter = new SAMFileWriterFactory().makeSAMWriter(alignedHeader, true, alignedSam)) {
            for (int i = 0; i < Math.max(firstOfPair.size(), secondOfPair.size()); ++i) {
                final HitSpec firstHitSpec = firstOfPair.isEmpty() ? null : firstOfPair.get(i);
                final HitSpec secondHitSpec = secondOfPair.isEmpty() ? null : secondOfPair.get(i);
                final SAMRecord first = makeRead(alignedHeader, firstUnmappedRec, firstHitSpec, true, i);
                final SAMRecord second = makeRead(alignedHeader, secondUnmappedRec, secondHitSpec, false, i);
                if (first != null && second != null) SamPairUtil.setMateInfo(first, second, alignedHeader);
                if (first != null) {
                    if (second == null) first.setMateUnmappedFlag(true);
                    alignedWriter.addAlignment(first);
                }
                if (second != null) {
                    if (first == null) second.setMateUnmappedFlag(true);
                    alignedWriter.addAlignment(second);
                }
            }
        }

        // Merge aligned file with original unmapped file.
        final File mergedSam = BaseTest.createTempFile("merged.", ".sam");

        doMergeAlignment(multiHitFilterUnmappedBam, Collections.singletonList(alignedSam),
                null, null, null, null,
                false, true, false, 1,
                "0", "1.0", "align!", "myAligner",
                fasta, mergedSam,
                null, null, null, null
        );

        SamAssertionUtils.assertSamValid(mergedSam);

        // Tally metrics and check for agreement with expected.
        final SamReader mergedReader = SamReaderFactory.makeDefault().open(mergedSam);
        int numFirst = 0;
        int numSecond = 0;
        Integer primaryHitIndex = null;
        int primaryMapq = 0;
        for (final SAMRecord rec : mergedReader) {
            if (rec.getFirstOfPairFlag()) ++numFirst;
            if (rec.getSecondOfPairFlag()) ++numSecond;
            if (!rec.getNotPrimaryAlignmentFlag()  && !rec.getReadUnmappedFlag()) {
                final Integer hitIndex = rec.getIntegerAttribute(SAMTag.HI.name());
                final int newHitIndex = (hitIndex == null? -1: hitIndex);
                if (primaryHitIndex == null) primaryHitIndex = newHitIndex;
                else Assert.assertEquals(newHitIndex, primaryHitIndex.intValue());
                primaryMapq += rec.getMappingQuality();
            }
        }
        Assert.assertEquals(numFirst, expectedNumFirst);
        Assert.assertEquals(numSecond, expectedNumSecond);
        Assert.assertEquals(primaryHitIndex, expectedPrimaryHitIndex);
        Assert.assertEquals(primaryMapq, expectedPrimaryMapq);
    }

    private SAMRecord makeRead(final SAMFileHeader alignedHeader, final SAMRecord unmappedRec, final HitSpec hitSpec,
                               final boolean firstOfPair, final int hitIndex) {
        if (hitSpec == null) return null;
        final SAMRecord rec = makeRead(alignedHeader, unmappedRec, hitSpec, hitIndex);
        rec.setReadPairedFlag(true);
        if (firstOfPair) {
            rec.setFirstOfPairFlag(true);
            rec.setAlignmentStart(hitIndex + 1);
        } else {
            rec.setSecondOfPairFlag(true);
            rec.setAlignmentStart(hitIndex + 201);
        }
        return rec;
    }

    private SAMRecord makeRead(final SAMFileHeader alignedHeader, final SAMRecord unmappedRec, final HitSpec hitSpec, final int hitIndex) {
        final SAMRecord rec = new SAMRecord(alignedHeader);
        rec.setReadName(unmappedRec.getReadName());
        rec.setReadBases(unmappedRec.getReadBases());
        rec.setBaseQualities(unmappedRec.getBaseQualities());
        rec.setMappingQuality(hitSpec.mapq);
        if (!hitSpec.primary) rec.setNotPrimaryAlignmentFlag(true);
        final Cigar cigar = new Cigar();
        final int readLength = rec.getReadLength();
        if (hitSpec.filtered) {
            // Add two insertions so alignment is filtered.
            cigar.add(new CigarElement(readLength-4, CigarOperator.M));
            cigar.add(new CigarElement(1, CigarOperator.I));
            cigar.add(new CigarElement(1, CigarOperator.M));
            cigar.add(new CigarElement(1, CigarOperator.I));
            cigar.add(new CigarElement(1, CigarOperator.M));
        } else {
            cigar.add(new CigarElement(readLength, CigarOperator.M));
        }
        rec.setCigar(cigar);

        rec.setReferenceName(bigSequenceName);
        rec.setAttribute(SAMTag.HI.name(), hitIndex);
        rec.setAlignmentStart(hitIndex + 1);

        return rec;
    }

    class HitSpec {
        final boolean primary;
        final boolean filtered;
        final int mapq;

        HitSpec(final boolean primary, final boolean filtered, final int mapq) {
            this.primary = primary;
            this.filtered = filtered;
            this.mapq = mapq;
        }
    }

    @DataProvider(name="testPairedMultiHitWithFilteringTestCases")
    public Object[][] testPairedMultiHitWithFilteringTestCases() {
        final ArrayList<Object[]> ret = new ArrayList<>();
        List<HitSpec> firstOfPair;
        List<HitSpec> secondOfPair;

        // One end aligned test cases...

        firstOfPair = new ArrayList<>();
        secondOfPair = Collections.emptyList();
        firstOfPair.add(new HitSpec(false, true, 10));
        firstOfPair.add(new HitSpec(true, true, 10));
        firstOfPair.add(new HitSpec(false, false, 9));
        ret.add(new Object[]{"Only 1st end has alignments, primary alignment is filtered, one 2ndary is filtered", firstOfPair, secondOfPair, -1, 1, 1, 9});

        secondOfPair = new ArrayList<>();
        firstOfPair = Collections.emptyList();
        secondOfPair.add(new HitSpec(false, false, 11));
        secondOfPair.add(new HitSpec(true, true, 10));
        secondOfPair.add(new HitSpec(false, true, 9));
        ret.add(new Object[]{"Only 2nd end has alignments, primary alignment is filtered, one 2ndary is filtered", firstOfPair, secondOfPair, -1, 1, 1, 11});

        firstOfPair = new ArrayList<>();
        secondOfPair = Collections.emptyList();
        firstOfPair.add(new HitSpec(false, true, 10));
        firstOfPair.add(new HitSpec(true, true, 10));
        firstOfPair.add(new HitSpec(false, true, 9));
        ret.add(new Object[]{"Only 1st end has alignments, all are filtered", firstOfPair, secondOfPair, null, 1, 1, 0});

        secondOfPair = new ArrayList<>();
        firstOfPair = Collections.emptyList();
        secondOfPair.add(new HitSpec(false, true, 10));
        secondOfPair.add(new HitSpec(true, true, 10));
        secondOfPair.add(new HitSpec(false, true, 9));
        ret.add(new Object[]{"Only 2nd end has alignments, all are filtered", firstOfPair, secondOfPair, null, 1, 1, 0});

        // Both ends aligned test cases...

        firstOfPair = new ArrayList<>();
        secondOfPair = new ArrayList<>();
        firstOfPair.add(new HitSpec(false, false, 10));
        firstOfPair.add(new HitSpec(true, true, 10));
        firstOfPair.add(new HitSpec(false, false, 9));
        secondOfPair.add(new HitSpec(false, false, 11));
        secondOfPair.add(new HitSpec(true, false, 30));
        secondOfPair.add(new HitSpec(false, false, 9));
        ret.add(new Object[]{"Both ends aligned, one end of primary filtered", firstOfPair, secondOfPair, -1, 3, 3, 30});

        firstOfPair = new ArrayList<>();
        secondOfPair = new ArrayList<>();
        firstOfPair.add(new HitSpec(false, false, 10));
        firstOfPair.add(new HitSpec(true, true, 10));
        firstOfPair.add(new HitSpec(false, true, 9));
        secondOfPair.add(new HitSpec(false, false, 11));
        secondOfPair.add(new HitSpec(true, false, 30));
        secondOfPair.add(new HitSpec(false, false, 9));
        ret.add(new Object[]{"Both ends aligned, one end of primary filtered, one end of secondary filtered", firstOfPair, secondOfPair, -1, 2, 3, 30});

        firstOfPair = new ArrayList<>();
        secondOfPair = new ArrayList<>();
        firstOfPair.add(new HitSpec(false, false, 9));
        firstOfPair.add(new HitSpec(true, true, 10));
        firstOfPair.add(new HitSpec(false, false, 10));
        secondOfPair.add(new HitSpec(false, false, 9));
        secondOfPair.add(new HitSpec(true, true, 30));
        secondOfPair.add(new HitSpec(false, false, 11));
        ret.add(new Object[]{"Both ends aligned, both ends of primary filtered", firstOfPair, secondOfPair, 1, 2, 2, 21});

        firstOfPair = new ArrayList<>();
        secondOfPair = new ArrayList<>();
        firstOfPair.add(new HitSpec(false, false, 9));
        firstOfPair.add(new HitSpec(true, true, 10));
        firstOfPair.add(new HitSpec(false, false, 10));
        secondOfPair.add(new HitSpec(false, false, 12));
        secondOfPair.add(new HitSpec(true, true, 30));
        secondOfPair.add(new HitSpec(false, false, 11));
        ret.add(new Object[]{"Both ends aligned, both ends of primary filtered, two secondary alignments with identical mapq",
                firstOfPair, secondOfPair, 1, 2, 2, 21});

        firstOfPair = new ArrayList<>();
        secondOfPair = new ArrayList<>();
        firstOfPair.add(new HitSpec(false, true, 9));
        firstOfPair.add(new HitSpec(true, true, 10));
        firstOfPair.add(new HitSpec(false, false, 10));
        secondOfPair.add(new HitSpec(false, false, 30));
        secondOfPair.add(new HitSpec(true, true, 30));
        secondOfPair.add(new HitSpec(false, false, 11));
        ret.add(new Object[]{"Both ends aligned, both ends of primary filtered, one end of secondary filtered, but with high mapq",
                firstOfPair, secondOfPair, -1, 2, 2, 30});

        firstOfPair = new ArrayList<>();
        secondOfPair = new ArrayList<>();
        firstOfPair.add(new HitSpec(false, true, 9));
        firstOfPair.add(new HitSpec(false, false, 10));
        firstOfPair.add(new HitSpec(true, false, 10));
        secondOfPair.add(new HitSpec(false, true, 30));
        secondOfPair.add(new HitSpec(false, false, 11));
        secondOfPair.add(new HitSpec(true, false, 30));
        ret.add(new Object[]{"Both ends aligned, secondary filtered on both ends", firstOfPair, secondOfPair, 1, 2, 2, 40});

        firstOfPair = new ArrayList<>();
        secondOfPair = new ArrayList<>();
        firstOfPair.add(new HitSpec(false, false, 9));
        firstOfPair.add(new HitSpec(false, false, 10));
        firstOfPair.add(new HitSpec(true, false, 10));
        secondOfPair.add(new HitSpec(false, true, 30));
        secondOfPair.add(new HitSpec(false, false, 11));
        secondOfPair.add(new HitSpec(true, false, 30));
        ret.add(new Object[]{"Both ends aligned, secondary filtered on second end", firstOfPair, secondOfPair, 2, 3, 2, 40});

        firstOfPair = new ArrayList<>();
        secondOfPair = new ArrayList<>();
        firstOfPair.add(new HitSpec(false, false, 9));
        firstOfPair.add(new HitSpec(false, true, 10));
        firstOfPair.add(new HitSpec(true, false, 10));
        secondOfPair.add(new HitSpec(false, false, 30));
        secondOfPair.add(new HitSpec(false, false, 11));
        secondOfPair.add(new HitSpec(true, false, 30));
        ret.add(new Object[]{"Both ends aligned, secondary filtered on first end", firstOfPair, secondOfPair, 2, 2, 3, 40});

        firstOfPair = new ArrayList<>();
        secondOfPair = new ArrayList<>();
        firstOfPair.add(new HitSpec(false, true, 9));
        firstOfPair.add(new HitSpec(false, true, 10));
        firstOfPair.add(new HitSpec(true, true, 10));
        secondOfPair.add(new HitSpec(false, true, 30));
        secondOfPair.add(new HitSpec(false, true, 11));
        secondOfPair.add(new HitSpec(true, true, 30));
        ret.add(new Object[]{"Both ends aligned, all filtered", firstOfPair, secondOfPair, null, 1, 1, 0});

        return ret.toArray(new Object[0][]);
    }




    /**
     * Read a single fragment read from a file, and create one or more aligned records for the read pair based on
     * the lists, merge with the original read, and assert expected results.
     * @param description
     * @param hitSpecs List that describes the aligned SAMRecords to create.
     * @param expectedPrimaryHitIndex Expected value for the HI tag in the primary alignment in the merged output.
     * @param expectedNumReads Expected number of SAMRecords in the merged output.
     * @param expectedPrimaryMapq Mapq of both ends of primary alignment in the merged output.
     * @throws Exception
     */
    @Test(dataProvider = "testFragmentMultiHitWithFilteringTestCases")
    public void testFragmentMultiHitWithFiltering(final String description, final List<HitSpec> hitSpecs,
                                                  final Integer expectedPrimaryHitIndex, final int expectedNumReads,
                                                  final int expectedPrimaryMapq) throws Exception {

        // Create the aligned file by copying bases, quals, readname from the unmapped read, and conforming to each HitSpec.
        try (final SAMRecordIterator unmappedSamFileIterator = SamReaderFactory.makeDefault().open(multiHitFilterFragmentUnmappedBam).iterator()) {
            final SAMRecord unmappedRec = unmappedSamFileIterator.next();
            final File alignedSam = BaseTest.createTempFile("aligned.", ".sam");
            final SAMFileHeader alignedHeader = new SAMFileHeader();
            alignedHeader.setSortOrder(SAMFileHeader.SortOrder.queryname);
            alignedHeader.setSequenceDictionary(SamReaderFactory.makeDefault().getFileHeader(sequenceDict).getSequenceDictionary());
            try (final SAMFileWriter alignedWriter = new SAMFileWriterFactory().makeSAMWriter(alignedHeader, true, alignedSam)) {
                for (int i = 0; i < hitSpecs.size(); ++i) {
                    final HitSpec hitSpec = hitSpecs.get(i);
                    final SAMRecord mappedRec = makeRead(alignedHeader, unmappedRec, hitSpec, i);
                    if (mappedRec != null) {
                        alignedWriter.addAlignment(mappedRec);
                    }
                }
            }
            // Merge aligned file with original unmapped file.
            final File mergedSam = BaseTest.createTempFile("merged.", ".sam");

            doMergeAlignment(multiHitFilterFragmentUnmappedBam, Collections.singletonList(alignedSam),
                    null, null, null, null,
                    false, true, false, 1,
                    "0", "1.0", "align!", "myAligner",
                    fasta, mergedSam,
                    null, null, null, null);

            SamAssertionUtils.assertSamValid(mergedSam);

            // Tally metrics and check for agreement with expected.
            final SamReader mergedReader = SamReaderFactory.makeDefault().open(mergedSam);
            int numReads = 0;
            Integer primaryHitIndex = null;
            int primaryMapq = 0;
            for (final SAMRecord rec : mergedReader) {
                ++numReads;
                if (!rec.getNotPrimaryAlignmentFlag()  && !rec.getReadUnmappedFlag()) {
                    final Integer hitIndex = rec.getIntegerAttribute(SAMTag.HI.name());
                    final int newHitIndex = (hitIndex == null? -1: hitIndex);
                    Assert.assertNull(primaryHitIndex);
                    primaryHitIndex = newHitIndex;
                    primaryMapq = rec.getMappingQuality();
                }
            }
            Assert.assertEquals(numReads, expectedNumReads);
            Assert.assertEquals(primaryHitIndex, expectedPrimaryHitIndex);
            Assert.assertEquals(primaryMapq, expectedPrimaryMapq);
        }

    }

    @DataProvider(name="testFragmentMultiHitWithFilteringTestCases")
    public Object[][] testFragmentMultiHitWithFilteringTestCases() {
        final ArrayList<Object[]> ret = new ArrayList<>();
        List<HitSpec> hitSpecs;

        hitSpecs = new ArrayList<>();
        hitSpecs.add(new HitSpec(false, true, 10));
        hitSpecs.add(new HitSpec(true, false, 8));
        hitSpecs.add(new HitSpec(false, false, 9));
        ret.add(new Object[]{"One secondary filtered", hitSpecs, -1, 2, 8});

        hitSpecs = new ArrayList<>();
        hitSpecs.add(new HitSpec(false, false, 10));
        hitSpecs.add(new HitSpec(true, true, 8));
        hitSpecs.add(new HitSpec(false, false, 11));
        ret.add(new Object[]{"Primary filtered", hitSpecs, -1, 2, 11});

        hitSpecs = new ArrayList<>();
        hitSpecs.add(new HitSpec(false, false, 11));
        hitSpecs.add(new HitSpec(true, true, 8));
        hitSpecs.add(new HitSpec(false, false, 11));
        ret.add(new Object[]{"Primary filtered, two secondaries with identical mapq", hitSpecs, -1, 2, 11});

        hitSpecs = new ArrayList<>();
        hitSpecs.add(new HitSpec(false, true, 10));
        hitSpecs.add(new HitSpec(true, true, 8));
        hitSpecs.add(new HitSpec(false, false, 9));
        ret.add(new Object[]{"Primary and one secondary filtered", hitSpecs, -1, 1, 9});

        hitSpecs = new ArrayList<>();
        hitSpecs.add(new HitSpec(false, true, 10));
        hitSpecs.add(new HitSpec(true, true, 8));
        hitSpecs.add(new HitSpec(false, true, 9));
        ret.add(new Object[]{"All filtered", hitSpecs, null, 1, 0});

        return ret.toArray(new Object[0][]);
    }

    /**
     * Confirm that paired reads are rejected by PrimaryAlignmentStrategy.EarliestFragment.
     */
    @Test(expectedExceptions = UnsupportedOperationException.class)
    public void testEarliestFragmentStrategyPaired() throws Exception {
        final File output = BaseTest.createTempFile("mergeTest", ".sam");

        final File unmappedSam = BaseTest.createTempFile("unmapped.", ".sam");
        final SAMFileWriterFactory factory = new SAMFileWriterFactory();
        final SAMFileHeader header = new SAMFileHeader();
        header.setSortOrder(SAMFileHeader.SortOrder.queryname);
        final String cigar = "16M";

        final SAMRecord firstOfPair = new SAMRecord(header);
        firstOfPair.setReadName("theRead");
        firstOfPair.setReadString("ACGTACGTACGTACGT");
        firstOfPair.setBaseQualityString("5555555555555555");
        firstOfPair.setReadUnmappedFlag(true);
        firstOfPair.setReadPairedFlag(true);
        firstOfPair.setFirstOfPairFlag(true);

        final SAMRecord secondOfPair = new SAMRecord(header);
        secondOfPair.setReadName("theRead");
        secondOfPair.setReadString("ACGTACGTACGTACGT");
        secondOfPair.setBaseQualityString("5555555555555555");
        secondOfPair.setReadUnmappedFlag(true);
        secondOfPair.setReadPairedFlag(true);
        secondOfPair.setSecondOfPairFlag(true);
        SamPairUtil.setMateInfo(firstOfPair, secondOfPair, header);

        try (final SAMFileWriter unmappedWriter = factory.makeSAMWriter(header, false, unmappedSam)) {
            unmappedWriter.addAlignment(firstOfPair);
            unmappedWriter.addAlignment(secondOfPair);
        }

        final File alignedSam = BaseTest.createTempFile("aligned.", ".sam");

        // Populate the header with SAMSequenceRecords
        header.getSequenceDictionary().addSequence(new SAMSequenceRecord("chr1", 1000000));

        // Create 2 alignments for each end of pair
        try (final SAMFileWriter alignedWriter = factory.makeSAMWriter(header, false, alignedSam)) {
            for (int i = 1; i <= 2; ++i) {
                final SAMRecord firstOfPairAligned = new SAMRecord(header);
                firstOfPairAligned.setReadName(firstOfPair.getReadName());
                firstOfPairAligned.setReadBases(firstOfPair.getReadBases());
                firstOfPairAligned.setBaseQualities(firstOfPair.getBaseQualities());
                firstOfPairAligned.setReferenceName("chr1");
                firstOfPairAligned.setAlignmentStart(i);
                firstOfPairAligned.setCigarString(cigar);
                firstOfPairAligned.setMappingQuality(100);
                firstOfPairAligned.setReadPairedFlag(true);
                firstOfPairAligned.setFirstOfPairFlag(true);
                firstOfPairAligned.setAttribute(SAMTag.HI.name(), i);

                final SAMRecord secondOfPairAligned = new SAMRecord(header);
                secondOfPairAligned.setReadName(secondOfPair.getReadName());
                secondOfPairAligned.setReadBases(secondOfPair.getReadBases());
                secondOfPairAligned.setBaseQualities(secondOfPair.getBaseQualities());
                secondOfPairAligned.setReferenceName("chr1");
                secondOfPairAligned.setAlignmentStart(i + 10);
                secondOfPairAligned.setCigarString(cigar);
                secondOfPairAligned.setMappingQuality(100);
                secondOfPairAligned.setReadPairedFlag(true);
                secondOfPairAligned.setSecondOfPairFlag(true);
                secondOfPairAligned.setAttribute(SAMTag.HI.name(), i);

                SamPairUtil.setMateInfo(firstOfPairAligned, secondOfPairAligned, header);

                alignedWriter.addAlignment(firstOfPairAligned);
                alignedWriter.addAlignment(secondOfPairAligned);
            }
        }

        doMergeAlignment(unmappedSam, Collections.singletonList(alignedSam),
                null, null, null, null,
                false, true, false, 1,
                "0", "1.0", "align!", "myAligner",
                fasta, output,
                SamPairUtil.PairOrientation.FR,
                MergeBamAlignment.PrimaryAlignmentStrategy.EarliestFragment,
                null, null);

        Assert.fail("Exception was not thrown");
    }

    /**
     * Various scenarios for EarliestFragmentStrategy.  Confirms that one of the expected ones is selected.
     * Note that there may be an arbitrary selection due to a tie.
     */
    @Test(dataProvider = "testEarliestFragmentStrategyDataProvider")
    public void testEarliestFragmentStrategy(final String testName, final MultipleAlignmentSpec[] specs) throws IOException {

        final File output = BaseTest.createTempFile(testName, ".sam");
        final File[] sams = createSamFilesToBeMerged(specs);

        doMergeAlignment(sams[0], Collections.singletonList(sams[1]),
                null, null, null, null,
                false, true, false, 1,
                "0", "1.0", "align!", "myAligner",
                fasta, output,
                SamPairUtil.PairOrientation.FR, MergeBamAlignment.PrimaryAlignmentStrategy.EarliestFragment,
                ONE_OF_THE_BEST_TAG,
                null
        );


        final SamReader mergedReader = SamReaderFactory.makeDefault().open(output);
        boolean seenPrimary = false;
        for (final SAMRecord rec : mergedReader) {
            if (!rec.getNotPrimaryAlignmentFlag()) {
                seenPrimary = true;
                final Integer oneOfTheBest = rec.getIntegerAttribute(ONE_OF_THE_BEST_TAG);
                Assert.assertEquals(oneOfTheBest, new Integer(1), "Read not marked as one of the best is primary: " + rec);
            }
        }
        CloserUtil.close(mergedReader);
        Assert.assertTrue(seenPrimary, "Never saw primary alignment");
    }

    @DataProvider(name = "testEarliestFragmentStrategyDataProvider")
    public Object[][] testEarliestFragmentStrategyDataProvider() {
        return new Object[][] {
                {"simpleForward", new MultipleAlignmentSpec[]{new MultipleAlignmentSpec("16M", false, 200, true)}},
                {"simpleReverse", new MultipleAlignmentSpec[]{new MultipleAlignmentSpec("16M", true, 200, true)}},
                {"2 forward one earlier", new MultipleAlignmentSpec[]{
                        new MultipleAlignmentSpec("1S15M", false, 200, false),
                        new MultipleAlignmentSpec("16M", false, 195, true)}},
                {"forward earlier than reverse", new MultipleAlignmentSpec[]{
                        new MultipleAlignmentSpec("16M", false, 200, true),
                        new MultipleAlignmentSpec("15M1S", true, 195, false)}},
                {"reverse earlier than forward", new MultipleAlignmentSpec[]{
                        new MultipleAlignmentSpec("2S14M", false, 200, false),
                        new MultipleAlignmentSpec("15M1S", true, 200, true)}},
                {"tie resolved via MAPQ", new MultipleAlignmentSpec[]{
                        new MultipleAlignmentSpec("1S15M", false, 200, false),
                        new MultipleAlignmentSpec("1S13M1S", false, 205, true),
                        new MultipleAlignmentSpec("15M1S", true, 200, false),
                        new MultipleAlignmentSpec("14M2S", true, 195, false)}},
                {"tie with same MAPQ resolved arbitrarily", new MultipleAlignmentSpec[]{
                        new MultipleAlignmentSpec("1S15M", false, 200, false),
                        new MultipleAlignmentSpec("1S13M1S", false, 205, true),
                        new MultipleAlignmentSpec("15M1S", true, 205, true),
                        new MultipleAlignmentSpec("14M2S", true, 195, false)}},
                {"one cigar with deletion and higher MAPQ", new MultipleAlignmentSpec[]{
                        new MultipleAlignmentSpec("16M", false, 200, false),
                        new MultipleAlignmentSpec("1D16M", false, 205, true)}},
                {"one cigar with deletion and lower MAPQ", new MultipleAlignmentSpec[]{
                        new MultipleAlignmentSpec("16M", false, 205, true),
                        new MultipleAlignmentSpec("1D16M", false, 200, false)}},
                {"Insertion makes alignment later", new MultipleAlignmentSpec[]{
                        new MultipleAlignmentSpec("16M", false, 200, true),
                        new MultipleAlignmentSpec("1I15M", false, 205, false)}},
                {"one cigar with deletion and higher MAPQ -- reverse", new MultipleAlignmentSpec[]{
                        new MultipleAlignmentSpec("16M", true, 200, false),
                        new MultipleAlignmentSpec("16M1D", true, 205, true)}},
                {"one cigar with deletion and lower MAPQ -- reverse", new MultipleAlignmentSpec[]{
                        new MultipleAlignmentSpec("16M", true, 205, true),
                        new MultipleAlignmentSpec("16M1D", true, 200, false)}},
                {"Insertion makes alignment later -- reverse", new MultipleAlignmentSpec[]{
                        new MultipleAlignmentSpec("16M", true, 200, true),
                        new MultipleAlignmentSpec("15M1I", true, 205, false)}}
        };
    }

    /**
     * @return a 2-element array in which the first element is the unmapped SAM, and the second the mapped SAM.
     */
    private File[] createSamFilesToBeMerged(final MultipleAlignmentSpec[] specs) {
        final File unmappedSam = BaseTest.createTempFile("unmapped.", ".sam");
        final SAMFileWriterFactory factory = new SAMFileWriterFactory();
        final SAMFileHeader header = new SAMFileHeader();
        header.setSortOrder(SAMFileHeader.SortOrder.queryname);
        final SAMRecord unmappedRecord = new SAMRecord(header);

        unmappedRecord.setReadName("theRead");
        unmappedRecord.setReadString("ACGTACGTACGTACGT");
        unmappedRecord.setBaseQualityString("5555555555555555");
        unmappedRecord.setReadUnmappedFlag(true);

        try (final SAMFileWriter unmappedWriter = factory.makeSAMWriter(header, false, unmappedSam)) {
            unmappedWriter.addAlignment(unmappedRecord);
        }

        final File alignedSam = BaseTest.createTempFile("aligned.", ".sam");

        final String sequence = "chr1";
        // Populate the header with SAMSequenceRecords
        header.getSequenceDictionary().addSequence(new SAMSequenceRecord(sequence, 1000000));

        try (final SAMFileWriter alignedWriter = factory.makeSAMWriter(header, false, alignedSam)){
            for (final MultipleAlignmentSpec spec : specs) {
                final SAMRecord alignedRecord = new SAMRecord(header);
                alignedRecord.setReadName(unmappedRecord.getReadName());
                alignedRecord.setReadBases(unmappedRecord.getReadBases());
                alignedRecord.setBaseQualities(unmappedRecord.getBaseQualities());
                alignedRecord.setReferenceName(sequence);
                alignedRecord.setAlignmentStart(1);
                alignedRecord.setReadNegativeStrandFlag(spec.reverseStrand);
                alignedRecord.setCigarString(spec.cigar);
                alignedRecord.setMappingQuality(spec.mapQ);
                if (spec.oneOfTheBest) {
                    alignedRecord.setAttribute(ONE_OF_THE_BEST_TAG, 1);
                }
                alignedWriter.addAlignment(alignedRecord);
            }
        }

        return new File[]{unmappedSam, alignedSam};
    }

    class MultipleAlignmentSpec {
        final String cigar;
        final boolean reverseStrand;
        final int mapQ;
        final boolean oneOfTheBest;

        MultipleAlignmentSpec(final String cigar, final boolean reverseStrand, final int mapQ, final boolean oneOfTheBest) {
            this.cigar = cigar;
            this.reverseStrand = reverseStrand;
            this.mapQ = mapQ;
            this.oneOfTheBest = oneOfTheBest;
        }
    }

    /**
     * Test that clipping of FR reads for fragments shorter than read length happens only when it should.
     */
    @Test
    public void testShortFragmentClipping() throws Exception {
        final File output = BaseTest.createTempFile("testShortFragmentClipping", ".sam");
        doMergeAlignment(clipTestUnmappedBam,
                Collections.singletonList(clipTestAlignedBam),
                null, null, null, null,
                false, true, false, 1,
                "0", "1.0", "align!", "myAligner",
                clippedReference, output,
                SamPairUtil.PairOrientation.FR, null,
                null, null);

        try (final SamReader result = SamReaderFactory.makeDefault().open(output)) {
            final Map<String, SAMRecord> firstReadEncountered = new HashMap<>();

            for (final SAMRecord rec : result) {
                final SAMRecord otherEnd = firstReadEncountered.get(rec.getReadName());
                if (otherEnd == null) {
                    firstReadEncountered.put(rec.getReadName(), rec);
                } else {
                    final int fragmentStart = Math.min(rec.getAlignmentStart(), otherEnd.getAlignmentStart());
                    final int fragmentEnd = Math.max(rec.getAlignmentEnd(), otherEnd.getAlignmentEnd());
                    final String[] readNameFields = rec.getReadName().split(":");
                    // Read name of each pair includes the expected fragment start and fragment end positions.
                    final int expectedFragmentStart = Integer.parseInt(readNameFields[1]);
                    final int expectedFragmentEnd = Integer.parseInt(readNameFields[2]);
                    Assert.assertEquals(fragmentStart, expectedFragmentStart, rec.getReadName());
                    Assert.assertEquals(fragmentEnd, expectedFragmentEnd, rec.getReadName());
                }
            }
        }
    }

    @Test(dataProvider="testBestFragmentMapqStrategy")
    public void testBestFragmentMapqStrategy(final String testName, final int[] firstMapQs, final int[] secondMapQs,
                                             final int expectedFirstMapq, final int expectedSecondMapq) throws Exception {
        testBestFragmentMapqStrategy(testName + "includeSecondary", firstMapQs, secondMapQs, true, expectedFirstMapq,
                expectedSecondMapq);
        testBestFragmentMapqStrategy(testName + "excludeSecondary", firstMapQs, secondMapQs, false, expectedFirstMapq,
                expectedSecondMapq);
    }

    private void testBestFragmentMapqStrategy(final String testName, final int[] firstMapQs, final int[] secondMapQs,
                                              final boolean includeSecondary, final int expectedFirstMapq,
                                              final int expectedSecondMapq) throws Exception {
        final File unmappedSam = BaseTest.createTempFile("unmapped.", ".sam");
        final SAMFileWriterFactory factory = new SAMFileWriterFactory();
        final SAMFileHeader header = new SAMFileHeader();
        header.setSortOrder(SAMFileHeader.SortOrder.queryname);

        final String readName = "theRead";
        final SAMRecord firstUnmappedRead = new SAMRecord(header);
        firstUnmappedRead.setReadName(readName);
        firstUnmappedRead.setReadString("ACGTACGTACGTACGT");
        firstUnmappedRead.setBaseQualityString("5555555555555555");
        firstUnmappedRead.setReadUnmappedFlag(true);
        firstUnmappedRead.setMateUnmappedFlag(true);
        firstUnmappedRead.setReadPairedFlag(true);
        firstUnmappedRead.setFirstOfPairFlag(true);

        final SAMRecord secondUnmappedRead = new SAMRecord(header);
        secondUnmappedRead.setReadName(readName);
        secondUnmappedRead.setReadString("TCGAACGTTCGAACTG");
        secondUnmappedRead.setBaseQualityString("6666666666666666");
        secondUnmappedRead.setReadUnmappedFlag(true);
        secondUnmappedRead.setMateUnmappedFlag(true);
        secondUnmappedRead.setReadPairedFlag(true);
        secondUnmappedRead.setSecondOfPairFlag(true);



        try (final SAMFileWriter unmappedWriter = factory.makeSAMWriter(header, false, unmappedSam)) {
            unmappedWriter.addAlignment(firstUnmappedRead);
            unmappedWriter.addAlignment(secondUnmappedRead);
        }

        final File alignedSam = BaseTest.createTempFile("aligned.", ".sam");

        final String sequence = "chr1";
        // Populate the header with SAMSequenceRecords
        header.getSequenceDictionary().addSequence(new SAMSequenceRecord(sequence, 1000000));

        try (final SAMFileWriter alignedWriter = factory.makeSAMWriter(header, false, alignedSam)) {
            addAlignmentsForBestFragmentMapqStrategy(alignedWriter, firstUnmappedRead, sequence, firstMapQs);
            addAlignmentsForBestFragmentMapqStrategy(alignedWriter, secondUnmappedRead, sequence, secondMapQs);
        }

        final File output = BaseTest.createTempFile("testBestFragmentMapqStrategy." + testName, ".sam");
        doMergeAlignment(unmappedSam, Collections.singletonList(alignedSam),
                null, null, null, null,
                false, true, false, 1,
                "0", "1.0", "align!", "myAligner",
                clippedReference, output,
                SamPairUtil.PairOrientation.FR,
                MergeBamAlignment.PrimaryAlignmentStrategy.BestEndMapq,
                null, includeSecondary
        );
        try (final SamReader reader = SamReaderFactory.makeDefault().open(output)) {

            int numFirstRecords = 0;
            int numSecondRecords = 0;
            int firstPrimaryMapq = -1;
            int secondPrimaryMapq = -1;
            for (final SAMRecord rec : reader) {
                Assert.assertTrue(rec.getReadPairedFlag());
                if (rec.getFirstOfPairFlag()) ++numFirstRecords;
                else if (rec.getSecondOfPairFlag()) ++numSecondRecords;
                else Assert.fail("unpossible!");
                if (!rec.getReadUnmappedFlag() && !rec.getNotPrimaryAlignmentFlag()) {
                    if (rec.getFirstOfPairFlag()) {
                        Assert.assertEquals(firstPrimaryMapq, -1);
                        firstPrimaryMapq = rec.getMappingQuality();
                    } else {
                        Assert.assertEquals(secondPrimaryMapq, -1);
                        secondPrimaryMapq = rec.getMappingQuality();
                    }
                } else if (rec.getNotPrimaryAlignmentFlag()) {
                    Assert.assertTrue(rec.getMateUnmappedFlag());
                }
            }
            Assert.assertEquals(firstPrimaryMapq, expectedFirstMapq);
            Assert.assertEquals(secondPrimaryMapq, expectedSecondMapq);
            if (!includeSecondary) {
                Assert.assertEquals(numFirstRecords, 1);
                Assert.assertEquals(numSecondRecords, 1);
            } else {
                // If no alignments for an end, there will be a single unmapped record
                Assert.assertEquals(numFirstRecords, Math.max(1, firstMapQs.length));
                Assert.assertEquals(numSecondRecords, Math.max(1, secondMapQs.length));
            }
        }
    }

    private void doMergeAlignment(final File unmappedBam, final List<File> alignedBams,
                                  final List<File> read1AlignedBams, final List<File> read2AlignedBams, final Integer read1Trim, final Integer read2Trim,
                                  final boolean alignReadsOnly, final boolean clipAdapters, final boolean isBisulfiteSequence, final int maxInsOrDels,
                                  final String progRecordId, final String progGroupVersion, final String progGroupCommandLine, final String progGroupName,
                                  final File refSeq, final File output,
                                  final SamPairUtil.PairOrientation expectedOrientation, final MergeBamAlignment.PrimaryAlignmentStrategy primaryAlignmentStrategy,
                                  final String attributesToRetain,
                                  final Boolean includeSecondary) {

        final List<String> args = new ArrayList<>(Arrays.asList(
                "--UNMAPPED_BAM", unmappedBam.getAbsolutePath(),
                "--ALIGNED_READS_ONLY", String.valueOf(alignReadsOnly),
                "--CLIP_ADAPTERS", String.valueOf(clipAdapters),
                "--IS_BISULFITE_SEQUENCE", String.valueOf(isBisulfiteSequence),
                "--MAX_INSERTIONS_OR_DELETIONS", String.valueOf(maxInsOrDels)));
        if (alignedBams != null) {
            for (final File alignedBam : alignedBams) {
                args.add("--ALIGNED_BAM"); args.add( alignedBam.getAbsolutePath());
            }
        }
        if (read1AlignedBams != null) {
            for (final File alignedBam : read1AlignedBams) {
                args.add("--READ1_ALIGNED_BAM"); args.add(  alignedBam.getAbsolutePath());
            }
        }
        if (read2AlignedBams != null) {
            for (final File alignedBam : read2AlignedBams) {
                args.add("--READ2_ALIGNED_BAM"); args.add(  alignedBam.getAbsolutePath());
            }
        }
        if (read1Trim != null) {
            args.add("--READ1_TRIM"); args.add(String.valueOf(read1Trim));
        }
        if (read2Trim != null) {
            args.add("--READ2_TRIM"); args.add(String.valueOf(read2Trim));
        }
        if (progRecordId != null) {
            args.add("--PROGRAM_RECORD_ID"); args.add( progRecordId);
        }
        if (progGroupVersion != null) {
            args.add("--PROGRAM_GROUP_VERSION"); args.add(progGroupVersion);
        }
        if (progGroupCommandLine != null) {
            args.add("--PROGRAM_GROUP_COMMAND_LINE"); args.add(progGroupCommandLine);
        }
        if (progGroupName != null) {
            args.add("--PROGRAM_GROUP_NAME");args.add(progGroupName);
        }
        args.add("--reference"); args.add(refSeq.getAbsolutePath());
        args.add("--output"); args.add( output.getAbsolutePath());

        if (expectedOrientation != null) {
            args.add("--EXPECTED_ORIENTATIONS"); args.add(expectedOrientation.toString());
        }
        if (primaryAlignmentStrategy != null) {
            args.add("--PRIMARY_ALIGNMENT_STRATEGY"); args.add( primaryAlignmentStrategy.toString());
        }
        if (attributesToRetain != null) {
            args.add("--ATTRIBUTES_TO_RETAIN"); args.add(attributesToRetain);
        }
        if (includeSecondary != null) {
            args.add("--INCLUDE_SECONDARY_ALIGNMENTS"); args.add(String.valueOf(includeSecondary));
        }
        runCommandLine(args);
    }

    private void addAlignmentsForBestFragmentMapqStrategy(
            final SAMFileWriter writer, final SAMRecord unmappedRecord, final String sequence, final int[] mapqs) {
        boolean reverse = false;
        int alignmentStart = 1;
        for (final int mapq : mapqs) {
            final SAMRecord alignedRecord = new SAMRecord(writer.getFileHeader());
            alignedRecord.setReadName(unmappedRecord.getReadName());
            alignedRecord.setReadBases(unmappedRecord.getReadBases());
            alignedRecord.setBaseQualities(unmappedRecord.getBaseQualities());
            alignedRecord.setReferenceName(sequence);
            alignedRecord.setAlignmentStart(alignmentStart);
            alignmentStart += 10; // Any old position will do
            alignedRecord.setReadNegativeStrandFlag(reverse);
            reverse = !reverse;
            alignedRecord.setCigarString(unmappedRecord.getReadBases().length + "M");
            alignedRecord.setMappingQuality(mapq);
            alignedRecord.setReadPairedFlag(unmappedRecord.getReadPairedFlag());
            alignedRecord.setFirstOfPairFlag(unmappedRecord.getFirstOfPairFlag());
            alignedRecord.setSecondOfPairFlag(unmappedRecord.getSecondOfPairFlag());
            alignedRecord.setMateUnmappedFlag(true);
            writer.addAlignment(alignedRecord);
        }
    }


    @DataProvider(name="testBestFragmentMapqStrategy")
    public Object[][] testBestFragmentMapqStrategyDataProvider() {
        /**
         *
         * @param testName
         * @param firstMapQs
         * @param secondMapQs
         * @param expectedFirstMapq
         * @param expectedSecondMapq
         */
        return new Object[][] {
                {"singleAlignmentFirstEnd", new int[]{12}, new int[0], 12, -1},
                {"singleAlignmentSecondEnd", new int[0], new int[]{12}, -1, 12},
                {"singleAlignmentBothEnd", new int[]{13}, new int[]{12}, 13, 12},
                {"multipleBothEnds1", new int[]{10, 10, 11, 11, 255, 0}, new int[]{14, 11, 1, 14}, 11, 14},
                {"multipleBothEnds2", new int[]{255, 0, 255}, new int[]{255, 255}, 255, 255},
                {"multipleFirstEnd", new int[]{10, 10, 11, 11, 12}, new int[0], 12, -1},
                {"multipleSecondEnd", new int[0], new int[]{10, 10, 0, 11, 12}, -1, 12},
                {"multipleFirstEndSingleSecondEnd", new int[]{10, 10, 11, 11, 12}, new int[]{255}, 12, 255},
                {"singleFirstEndMultipleSecondEnd", new int[]{0}, new int[]{10, 10, 11, 11, 12}, 0, 12},
        };
    }

    @Test(dataProvider = "testMostDistantStrategy")
    public void testMostDistantStrategy(final String testName,
                                        final MostDistantStrategyAlignmentSpec[] firstEndSpecs,
                                        final MostDistantStrategyAlignmentSpec[] secondEndSpecs) throws Exception {
        testMostDistantStrategy(testName +".includeSecondary", true, firstEndSpecs, secondEndSpecs);
        testMostDistantStrategy(testName +".excludeSecondary", false, firstEndSpecs, secondEndSpecs);
    }

    public void testMostDistantStrategy(final String testName, final boolean includeSecondary,
                                        final MostDistantStrategyAlignmentSpec[] firstEndSpecs,
                                        final MostDistantStrategyAlignmentSpec[] secondEndSpecs) throws Exception {

        final File unmappedSam = BaseTest.createTempFile("unmapped.", ".sam");
        final SAMFileWriterFactory factory = new SAMFileWriterFactory();
        final SAMFileHeader header = new SAMFileHeader();
        header.setSortOrder(SAMFileHeader.SortOrder.queryname);

        final String readName = "theRead";
        final SAMRecord firstUnmappedRead = new SAMRecord(header);
        firstUnmappedRead.setReadName(readName);
        firstUnmappedRead.setReadString("ACGT");
        firstUnmappedRead.setBaseQualityString("5555");
        firstUnmappedRead.setReadUnmappedFlag(true);
        firstUnmappedRead.setMateUnmappedFlag(true);
        firstUnmappedRead.setReadPairedFlag(true);
        firstUnmappedRead.setFirstOfPairFlag(true);

        final SAMRecord secondUnmappedRead = new SAMRecord(header);
        secondUnmappedRead.setReadName(readName);
        secondUnmappedRead.setReadString("TCGA");
        secondUnmappedRead.setBaseQualityString("6666");
        secondUnmappedRead.setReadUnmappedFlag(true);
        secondUnmappedRead.setMateUnmappedFlag(true);
        secondUnmappedRead.setReadPairedFlag(true);
        secondUnmappedRead.setSecondOfPairFlag(true);

        try (final SAMFileWriter unmappedWriter = factory.makeSAMWriter(header, false, unmappedSam)) {
            unmappedWriter.addAlignment(firstUnmappedRead);
            unmappedWriter.addAlignment(secondUnmappedRead);
        }

        final File alignedSam = BaseTest.createTempFile("aligned.", ".sam");

        header.setSequenceDictionary(SamReaderFactory.makeDefault().getFileHeader(sequenceDict).getSequenceDictionary());

        String expectedFirstPrimarySequence = null;
        int expectedFirstPrimaryAlignmentStart = -1;
        String expectedSecondPrimarySequence = null;
        int expectedSecondPrimaryAlignmentStart = -1;

        try (final SAMFileWriter alignedWriter = factory.makeSAMWriter(header, false, alignedSam)) {

            // Semi-randomly make the reads align to forward or reverse strand.
            boolean reverse = false;
            for (final MostDistantStrategyAlignmentSpec spec : firstEndSpecs) {
                addAlignmentForMostStrategy(alignedWriter, firstUnmappedRead, spec, reverse);
                reverse = !reverse;
                if (spec.expectedPrimary) {
                    expectedFirstPrimarySequence = spec.sequence;
                    expectedFirstPrimaryAlignmentStart = spec.alignmentStart;
                }
            }
            for (final MostDistantStrategyAlignmentSpec spec : secondEndSpecs) {
                addAlignmentForMostStrategy(alignedWriter, secondUnmappedRead, spec, reverse);
                reverse = !reverse;
                if (spec.expectedPrimary) {
                    expectedSecondPrimarySequence = spec.sequence;
                    expectedSecondPrimaryAlignmentStart = spec.alignmentStart;
                }
            }
        }
        final File output = BaseTest.createTempFile("testMostDistantStrategy." + testName, ".sam");
        doMergeAlignment(unmappedSam, Collections.singletonList(alignedSam),
                null, null, null, null,
                false, true, false, 1,
                "0", "1.0", "align!", "myAligner",
                fasta, output,
                SamPairUtil.PairOrientation.FR, MergeBamAlignment.PrimaryAlignmentStrategy.MostDistant,
                null, includeSecondary);

        final SamReader reader = SamReaderFactory.makeDefault().open(output);
        int numFirstRecords = 0;
        int numSecondRecords = 0;
        String firstPrimarySequence = null;
        int firstPrimaryAlignmentStart = -1;
        String secondPrimarySequence = null;
        int secondPrimaryAlignmentStart = -1;
        for (final SAMRecord rec: reader) {
            Assert.assertTrue(rec.getReadPairedFlag());
            if (rec.getFirstOfPairFlag()) ++numFirstRecords;
            else if (rec.getSecondOfPairFlag()) ++ numSecondRecords;
            else Assert.fail("unpossible!");
            if (!rec.getReadUnmappedFlag() && !rec.getNotPrimaryAlignmentFlag()) {
                if (rec.getFirstOfPairFlag()) {
                    Assert.assertEquals(firstPrimaryAlignmentStart, -1);
                    firstPrimarySequence = rec.getReferenceName();
                    firstPrimaryAlignmentStart = rec.getAlignmentStart();
                } else {
                    Assert.assertEquals(secondPrimaryAlignmentStart, -1);
                    secondPrimarySequence = rec.getReferenceName();
                    secondPrimaryAlignmentStart = rec.getAlignmentStart();
                }
            } else if (rec.getNotPrimaryAlignmentFlag()) {
                Assert.assertTrue(rec.getMateUnmappedFlag());
            }
        }
        CloserUtil.close(reader);
        Assert.assertEquals(firstPrimarySequence, expectedFirstPrimarySequence);
        Assert.assertEquals(firstPrimaryAlignmentStart, expectedFirstPrimaryAlignmentStart);
        Assert.assertEquals(secondPrimarySequence, expectedSecondPrimarySequence);
        Assert.assertEquals(secondPrimaryAlignmentStart, expectedSecondPrimaryAlignmentStart);
        if (!includeSecondary) {
            Assert.assertEquals(numFirstRecords, 1);
            Assert.assertEquals(numSecondRecords, 1);
        } else {
            // If no alignments for an end, there will be a single unmapped record
            Assert.assertEquals(numFirstRecords, Math.max(1, firstEndSpecs.length));
            Assert.assertEquals(numSecondRecords, Math.max(1, secondEndSpecs.length));
        }
    }

    private void addAlignmentForMostStrategy(
            final SAMFileWriter writer, final SAMRecord unmappedRecord, final MostDistantStrategyAlignmentSpec spec,
            final boolean reverse) {
        final SAMRecord alignedRecord = new SAMRecord(writer.getFileHeader());
        alignedRecord.setReadName(unmappedRecord.getReadName());
        alignedRecord.setReadBases(unmappedRecord.getReadBases());
        alignedRecord.setBaseQualities(unmappedRecord.getBaseQualities());
        alignedRecord.setReferenceName(spec.sequence);
        alignedRecord.setAlignmentStart(spec.alignmentStart);
        alignedRecord.setReadNegativeStrandFlag(reverse);
        alignedRecord.setCigarString(unmappedRecord.getReadBases().length + "M");
        alignedRecord.setMappingQuality(spec.mapQ);
        alignedRecord.setReadPairedFlag(unmappedRecord.getReadPairedFlag());
        alignedRecord.setFirstOfPairFlag(unmappedRecord.getFirstOfPairFlag());
        alignedRecord.setSecondOfPairFlag(unmappedRecord.getSecondOfPairFlag());
        alignedRecord.setMateUnmappedFlag(true);
        writer.addAlignment(alignedRecord);
    }

    private static class MostDistantStrategyAlignmentSpec {
        public final boolean expectedPrimary;
        public final String sequence;
        public final int alignmentStart;
        public final int mapQ;

        private MostDistantStrategyAlignmentSpec(final boolean expectedPrimary, final String sequence,
                                                 final int alignmentStart, final int mapQ) {
            this.expectedPrimary = expectedPrimary;
            this.sequence = sequence;
            this.alignmentStart = alignmentStart;
            this.mapQ = mapQ;
        }

    }

    @DataProvider(name="testMostDistantStrategy")
    public Object[][] testMostDistantStrategyDataProvider() {
        /**
         * @param testName
         * @param firstEndSpecs
         * @param secondEndSpecs
         */
        return new Object[][] {
                {
                        // There are two ties: {chr1:1 - chr1:89} and {chr4:2 -chr4:90}
                        // That are disambiguated by MAPQ.
                        "multipleAlignmentsBothEnds",
                        new MostDistantStrategyAlignmentSpec[] {
                                new MostDistantStrategyAlignmentSpec(false, "chr1", 1, 20),
                                new MostDistantStrategyAlignmentSpec(false, "chr2", 1, 20),
                                new MostDistantStrategyAlignmentSpec(false, "chr3", 1, 20),
                                new MostDistantStrategyAlignmentSpec(true, "chr4", 2, 25),
                                new MostDistantStrategyAlignmentSpec(false, "chr1", 89, 20),
                                new MostDistantStrategyAlignmentSpec(false, "chr2", 2, 20),
                                new MostDistantStrategyAlignmentSpec(false, "chr3", 3, 20),
                                new MostDistantStrategyAlignmentSpec(false, "chr4", 4, 20),
                        },
                        new MostDistantStrategyAlignmentSpec[] {
                                new MostDistantStrategyAlignmentSpec(false, "chr1", 1, 20),
                                new MostDistantStrategyAlignmentSpec(false, "chr2", 2, 20),
                                new MostDistantStrategyAlignmentSpec(false, "chr3", 3, 20),
                                new MostDistantStrategyAlignmentSpec(true, "chr4", 90, 19),
                                new MostDistantStrategyAlignmentSpec(false, "chr1", 5, 20),
                                new MostDistantStrategyAlignmentSpec(false, "chr2", 6, 20),
                        },
                },
                {
                        "multipleAlignmentsAllChimeric",
                        new MostDistantStrategyAlignmentSpec[] {
                                new MostDistantStrategyAlignmentSpec(false, "chr1", 1, 20),
                                new MostDistantStrategyAlignmentSpec(false, "chr2", 2, 20),
                                new MostDistantStrategyAlignmentSpec(false, "chr1", 3, 20),
                                new MostDistantStrategyAlignmentSpec(true, "chr2", 4, 21),
                                new MostDistantStrategyAlignmentSpec(false, "chr1", 10, 20),
                                new MostDistantStrategyAlignmentSpec(false, "chr2", 50, 20),
                                new MostDistantStrategyAlignmentSpec(false, "chr1", 60, 20),
                                new MostDistantStrategyAlignmentSpec(false, "chr2", 70, 20),
                        },
                        new MostDistantStrategyAlignmentSpec[] {
                                new MostDistantStrategyAlignmentSpec(false, "chr3", 1, 20),
                                new MostDistantStrategyAlignmentSpec(false, "chr4", 2, 20),
                                new MostDistantStrategyAlignmentSpec(false, "chr3", 3, 20),
                                new MostDistantStrategyAlignmentSpec(true, "chr4", 11, 25),
                                new MostDistantStrategyAlignmentSpec(false, "chr3", 50, 20),
                                new MostDistantStrategyAlignmentSpec(false, "chr4", 60, 20),
                        },
                },
                {
                        "multipleAlignmentsFirstEnd",
                        new MostDistantStrategyAlignmentSpec[] {
                                new MostDistantStrategyAlignmentSpec(false, "chr1", 10, 20),
                                new MostDistantStrategyAlignmentSpec(false, "chr2", 10, 20),
                                new MostDistantStrategyAlignmentSpec(false, "chr3", 10, 20),
                                new MostDistantStrategyAlignmentSpec(true, "chr4", 20, 25),
                                new MostDistantStrategyAlignmentSpec(false, "chr1", 80, 20),
                                new MostDistantStrategyAlignmentSpec(false, "chr2", 20, 20),
                                new MostDistantStrategyAlignmentSpec(false, "chr3", 30, 20),
                                new MostDistantStrategyAlignmentSpec(false, "chr4", 40, 20),
                        },
                        new MostDistantStrategyAlignmentSpec[0]
                },
                {
                        "multipleAlignmentsSecondEnd",
                        new MostDistantStrategyAlignmentSpec[0],
                        new MostDistantStrategyAlignmentSpec[] {
                                new MostDistantStrategyAlignmentSpec(false, "chr1", 10, 20),
                                new MostDistantStrategyAlignmentSpec(false, "chr2", 10, 20),
                                new MostDistantStrategyAlignmentSpec(false, "chr3", 10, 20),
                                new MostDistantStrategyAlignmentSpec(true, "chr4", 20, 25),
                                new MostDistantStrategyAlignmentSpec(false, "chr1", 80, 20),
                                new MostDistantStrategyAlignmentSpec(false, "chr2", 20, 20),
                                new MostDistantStrategyAlignmentSpec(false, "chr3", 30, 20),
                                new MostDistantStrategyAlignmentSpec(false, "chr4", 40, 20),
                        },
                },
        };
    }

    @Test
    public void testHelp() {
        ByteArrayOutputStream baos = new ByteArrayOutputStream();
        PrintStream ps= new PrintStream(baos);
        System.setErr(ps);
        runCommandLine(new String[]{
                "--help",
        });
        final String s = baos.toString();
        String[] lines = s.split("\n");
        for (int i = 0; i < lines.length; i++) {
            Assert.assertFalse(lines[i].contains("READ2_ALIGNED_BAM (R2_ALIGNED) READ2_ALIGNED_BAM"), lines[i]); //Same option twice!
        }
    }
}
