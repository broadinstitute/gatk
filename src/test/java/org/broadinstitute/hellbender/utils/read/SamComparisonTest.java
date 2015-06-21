package org.broadinstitute.hellbender.utils.read;

import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

public final class SamComparisonTest {
    private static final File TEST_FILES_DIR = new File("src/test/resources/org/broadinstitute/hellbender/tools/picard/sam/CompareSAMs");

    private void testHelper(final String f1, final String f2,
                            final int expectedMatch, final int expectedDiffer, final int expectedUnmappedBoth,
                            final int expectedUnmappedLeft, final int expectedUnmappedRight, final int expectedMissingLeft,
                            final int expectedMissingRight, final boolean areEqual) throws IOException {

        SamReaderFactory factory = SamReaderFactory.makeDefault();
        File sam1 = new File(TEST_FILES_DIR, f1);
        File sam2 = new File(TEST_FILES_DIR, f2);

        try (final SamReader reader1 = factory.open(sam1);
             final SamReader reader2 = factory.open(sam2)) {
            final SamComparison comparison = new SamComparison(reader1, reader2);
            Assert.assertEquals(areEqual, comparison.areEqual());
            Assert.assertEquals(expectedMatch, comparison.getMappingsMatch());
            Assert.assertEquals(expectedDiffer, comparison.getMappingsDiffer());
            Assert.assertEquals(expectedUnmappedBoth, comparison.getUnmappedBoth());
            Assert.assertEquals(expectedUnmappedLeft, comparison.getUnmappedLeft());
            Assert.assertEquals(expectedUnmappedRight, comparison.getUnmappedRight());
            Assert.assertEquals(expectedMissingLeft, comparison.getMissingLeft());
            Assert.assertEquals(expectedMissingRight, comparison.getMissingRight());
        }

        // now reverse the comparison
        try (final SamReader reader1 = factory.open(sam1);
             final SamReader reader2 = factory.open(sam2)) {
            final SamComparison comparison = new SamComparison(reader2, reader1);
            Assert.assertEquals(areEqual, comparison.areEqual());
            Assert.assertEquals(expectedMatch, comparison.getMappingsMatch());
            Assert.assertEquals(expectedDiffer, comparison.getMappingsDiffer());
            Assert.assertEquals(expectedUnmappedBoth, comparison.getUnmappedBoth());
            Assert.assertEquals(expectedUnmappedRight, comparison.getUnmappedLeft());
            Assert.assertEquals(expectedUnmappedLeft, comparison.getUnmappedRight());
            Assert.assertEquals(expectedMissingRight, comparison.getMissingLeft());
            Assert.assertEquals(expectedMissingLeft, comparison.getMissingRight());
        }
    }

    @Test
    public void testSortsDifferent() throws IOException {
        // should be all 0's because program should return before comparing any alignments
        testHelper("genomic_sorted.sam", "unsorted.sam", 0, 0, 0, 0, 0, 0, 0, false);
    }

    @Test
    public void testSequenceDictionariesDifferent1() throws IOException {
        testHelper("genomic_sorted.sam", "chr21.sam", 0, 0, 0, 0, 0, 0, 0, false);
    }

    @Test
    public void testSequenceDictionariesDifferent2() throws IOException {
        testHelper("genomic_sorted.sam", "bigger_seq_dict.sam", 0, 0, 0, 0, 0, 0, 0, false);
    }

    @Test
    public void testBiggerSequenceDictionaries() throws IOException {
        testHelper("bigger_seq_dict.sam", "bigger_seq_dict.sam", 2, 0, 0, 0, 0, 0, 0, true);
    }

    @Test
    public void testIdentical() throws IOException {
        testHelper("genomic_sorted.sam", "genomic_sorted.sam", 2, 0, 0, 0, 0, 0, 0, true);
    }

    @Test
    public void testHasNonPrimary() throws IOException {
        testHelper("genomic_sorted.sam", "has_non_primary.sam", 2, 0, 0, 0, 0, 0, 0, true);
    }

    @Test
    public void testMoreOnOneSide() throws IOException {
        testHelper("genomic_sorted_5.sam", "genomic_sorted_5_plus.sam", 3, 2, 0, 0, 0, 3, 0, false);
    }

    @Test
    public void testGroupWithSameCoordinate() throws IOException {
        testHelper("group_same_coord.sam", "group_same_coord_diff_order.sam", 3, 0, 0, 0, 0, 1, 2, false);
    }

    @Test
    public void testGroupWithSameCoordinateAndNoMatchInOther() throws IOException {
        testHelper("group_same_coord.sam", "diff_coords.sam", 0, 5, 0, 0, 0, 0, 0, false);
    }

    @Test
    public void testUnmapped1() throws IOException {
        testHelper("genomic_sorted.sam", "unmapped_first.sam", 1, 0, 0, 0, 1, 0, 0, false);
    }

    @Test
    public void testUnmapped2() throws IOException {
        testHelper("genomic_sorted.sam", "unmapped_second.sam", 1, 0, 0, 0, 1, 0, 0, false);
    }

    @Test
    public void testUnmapped3() throws IOException {
        testHelper("unmapped_first.sam", "unmapped_second.sam", 0, 0, 0, 1, 1, 0, 0, false);
    }

    @Test
    public void testUnmapped4() throws IOException {
        testHelper("unmapped_first.sam", "unmapped_first.sam", 1, 0, 1, 0, 0, 0, 0, true);
    }

    @Test
    public void testDuplicateReadNameAndCoordinate() throws IOException {
        testHelper("same_qname_and_pos.sam", "same_qname_and_pos.sam", 1, 0, 1, 0, 0, 0, 0, true);
    }

    @Test
    public void testQuerynameSomeDifferentCoords() throws IOException {
        testHelper("qname_sorted_v1.sam", "qname_sorted_v2.sam", 2, 2, 0, 0, 0, 0, 0, false);
    }

    @Test
    public void testQuerynameMissingFirstOfPair() throws IOException {
        testHelper("qname_sorted_v1.sam", "qname_sorted_v1_second_of_pair_only.sam", 2, 0, 0, 0, 0, 0, 2, false);
    }

    @Test
    public void testUnsorted() throws IOException {
        testHelper("unsorted.sam", "unsorted_bigger.sam", 2, 0, 0, 0, 0, 1, 0, false);
    }
}
