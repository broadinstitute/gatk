package org.broadinstitute.hellbender.engine.filters;

import htsjdk.samtools.*;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.clipping.ReadClipperTestUtils;
import org.broadinstitute.hellbender.utils.read.ArtificialSAMUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.lang.reflect.Method;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

import static org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary.*;

/**
 * Tests for the read filter library.
 */
public final class ReadFilterLibraryUnitTest {
    private static final int CHR_COUNT = 1;
    private static final int CHR_START = 1;
    private static final int CHR_SIZE = 1000;
    private static final int GROUP_COUNT = 5;

    private final SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeaderWithGroups(CHR_COUNT, CHR_START, CHR_SIZE, GROUP_COUNT);

    /**
     * Creates a read record.
     *
     * @param cigar the new record CIGAR.
     * @param group the new record group index that must be in the range \
     *              [0,{@link #GROUP_COUNT})
     * @param reference the reference sequence index (0-based)
     * @param start the start position of the read alignment in the reference
     *              (1-based)
     * @return never <code>null</code>
     */
    private SAMRecord createRead(final Cigar cigar, final int group, final int reference, final int start) {
        final SAMRecord record = ArtificialSAMUtils.createArtificialRead(cigar);
        record.setHeader(header);
        record.setAlignmentStart(start);
        record.setReferenceIndex(reference);
        record.setAttribute(SAMTag.RG.toString(), header.getReadGroups().get(group).getReadGroupId());
        return record;

    }

    @Test
    public void testCheckSeqStored() {
        final SAMRecord goodRead = ArtificialSAMUtils.createArtificialRead(new byte[]{(byte) 'A'}, new byte[]{(byte) 'A'}, "1M");
        final SAMRecord badRead = ArtificialSAMUtils.createArtificialRead(new byte[]{}, new byte[]{}, "1M");
        badRead.setReadString("*");

        Assert.assertTrue(SEQ_IS_STORED.test(goodRead));
        Assert.assertFalse(SEQ_IS_STORED.test(badRead));
    }

    @Test(dataProvider = "UnsupportedCigarOperatorDataProvider")
    public void testCigarNOperatorFilter(String cigarString) {

        final ReadFilter filter = ReadFilterLibrary.WELLFORMED;
        final SAMRecord read = buildSAMRecord(cigarString);
        final boolean containsN = cigarString.contains("N");
        Assert.assertEquals(containsN, !filter.test(read), cigarString);
    }

    private SAMRecord buildSAMRecord(final String cigarString) {
        final Cigar nContainingCigar = TextCigarCodec.decode(cigarString);
        return this.createRead(nContainingCigar, 1, 0, 10);
    }


    @DataProvider(name = "UnsupportedCigarOperatorDataProvider")
    public Iterator<Object[]> unsupportedOperatorDataProvider(final Method testMethod) {
        /**
         * Cigar test data for unsupported operator test.
         * Each element of this array corresponds to a test case. In turn the first element of the test case array is the
         * Cigar string for that test case and the second indicates whether it should be filtered due to the presence of a
         * unsupported operator
         */
        final String[] TEST_CIGARS = {
                "101M10D20I10M",
                "6M14N5M",
                "1N",
                "101M",
                "110N",
                "2N4M",
                "4M2N",
                "3M1I1M",
                "1M2I2M",
                "1M10N1I1M",
                "1M1I1D",
                "11N12M1I34M12N"
        };
        final List<Object[]> result = new LinkedList<>();
        for (final String cigarString : TEST_CIGARS) {
            result.add(new Object[]{cigarString});
        }
        return result.iterator();
    }

    private SAMRecord simpleGoodRead() {
        String cigarString = "101M";
        final Cigar nContainingCigar = TextCigarCodec.decode(cigarString);
        SAMRecord read =  createRead(nContainingCigar, 1, 0, 10);
        read.setMappingQuality(50);
        return read;
    }

    @Test
    public void passesAllFilters() {
        SAMRecord read = simpleGoodRead();
        Assert.assertTrue(MAPPED.test(read), "MAPPED " + read.toString());
        Assert.assertTrue(PRIMARY_ALIGNMENT.test(read), "PRIMARY_ALIGNMENT " + read.toString());
        Assert.assertTrue(NOT_DUPLICATE.test(read), "NOT_DUPLICATE " + read.toString());
        Assert.assertTrue(PASSES_VENDOR_QUALITY_CHECK.test(read), "PASSES_VENDOR_QUALITY_CHECK " + read.toString());
        Assert.assertTrue(MAPPING_QUALITY_AVAILABLE.test(read), "MAPPING_QUALITY_AVAILABLE " + read.toString());
        Assert.assertTrue(MAPPING_QUALITY_NOT_ZERO.test(read), "MAPPING_QUALITY_NOT_ZERO " + read.toString());
        Assert.assertTrue(VALID_ALIGNMENT_START.test(read), "VALID_ALIGNMENT_START " + read.toString());
        Assert.assertTrue(VALID_ALIGNMENT_END.test(read), "VALID_ALIGNMENT_END " + read.toString());
        Assert.assertTrue(ALIGNMENT_AGREES_WITH_HEADER.test(read), "ALIGNMENT_AGREES_WITH_HEADER " + read.toString());
        Assert.assertTrue(HAS_READ_GROUP.test(read), "HAS_READ_GROUP " + read.toString());
        Assert.assertTrue(HAS_MATCHING_BASES_AND_QUALS.test(read), "HAS_MATCHING_BASES_AND_QUALS " + read.toString());
        Assert.assertTrue(SEQ_IS_STORED.test(read), "SEQ_IS_STORED " + read.toString());
        Assert.assertTrue(CIGAR_IS_SUPPORTED.test(read), "CIGAR_IS_SUPPORTED " + read.toString());

        Assert.assertTrue(WELLFORMED.test(read), "WELLFORMED " + read.toString());
    }

    @Test
    public void failsMAPPED_flag() {
        SAMRecord read = simpleGoodRead();
        read.setReadUnmappedFlag(true);
        Assert.assertFalse(MAPPED.test(read), read.toString());
    }

    @Test
    public void failsMAPPED_alignmenStart() {
        SAMRecord read = simpleGoodRead();
        read.setAlignmentStart(SAMRecord.NO_ALIGNMENT_START);
        Assert.assertFalse(MAPPED.test(read), read.toString());
    }

    @Test
    public void failsNOT_DUPLICATE() {
        SAMRecord read = simpleGoodRead();
        read.setDuplicateReadFlag(true);
        Assert.assertFalse(NOT_DUPLICATE.test(read), read.toString());
    }

    @Test
    public void failsPASSES_VENDOR_QUALITY_CHECK() {
        SAMRecord read = simpleGoodRead();
        read.setReadFailsVendorQualityCheckFlag(true);
        Assert.assertFalse(PASSES_VENDOR_QUALITY_CHECK.test(read), read.toString());
    }

    @Test
    public void failsMAPPING_QUALITY_AVAILABLE() {
        SAMRecord read = simpleGoodRead();
        read.setMappingQuality(QualityUtils.MAPPING_QUALITY_UNAVAILABLE);
        Assert.assertFalse(MAPPING_QUALITY_AVAILABLE.test(read), read.toString());
    }

    @Test
    public void failsMAPPING_QUALITY_NOT_ZERO() {
        SAMRecord read = simpleGoodRead();
        read.setMappingQuality(0);
        Assert.assertFalse(MAPPING_QUALITY_NOT_ZERO.test(read), read.toString());
    }

    @Test
    public void failsVALID_ALIGNMENT_START_case1() {
        SAMRecord read = simpleGoodRead();
        read.setAlignmentStart(SAMRecord.NO_ALIGNMENT_START);
        Assert.assertFalse(VALID_ALIGNMENT_START.test(read), read.toString());
    }

    @Test
    public void failsVALID_ALIGNMENT_START_case2() {
        SAMRecord read = simpleGoodRead();
        read.setAlignmentStart(-1);
        Assert.assertFalse(VALID_ALIGNMENT_START.test(read), read.toString());
    }

    @Test
    public void failsALIGNMENT_AGREES_WITH_HEADER_case1() {
        SAMRecord read = simpleGoodRead();
        read.setReferenceIndex(SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX);
        read.setAlignmentStart(10);
        Assert.assertFalse(ALIGNMENT_AGREES_WITH_HEADER.test(read), read.toString());
    }

    @Test
    public void failsALIGNMENT_AGREES_WITH_HEADER_case2() {
        SAMRecord read = simpleGoodRead();
        final int length = read.getHeader().getSequence(0).getSequenceLength();
        read.setAlignmentStart(length + 10);
        Assert.assertFalse(ALIGNMENT_AGREES_WITH_HEADER.test(read), read.toString());
    }

    @Test
    public void failsHAS_READ_GROUP() {
        SAMRecord read = simpleGoodRead();
        read.setAttribute(SAMTag.RG.name(), null);
        Assert.assertFalse(HAS_READ_GROUP.test(read), read.toString());
    }

    @Test
    public void failsHAS_MATCHING_BASES_AND_QUALS() {
        SAMRecord read = simpleGoodRead();
        read.setBaseQualities(new byte[]{1, 2, 3});
        Assert.assertFalse(HAS_MATCHING_BASES_AND_QUALS.test(read), read.toString());
    }

    @Test
    public void failsSEQ_IS_STORED() {
        SAMRecord read = simpleGoodRead();
        read.setReadBases(SAMRecord.NULL_SEQUENCE);
        Assert.assertFalse(SEQ_IS_STORED.test(read), read.toString());
    }

    @Test
    public void failsCIGAR_IS_SUPPORTED() {
        SAMRecord read = simpleGoodRead();
        read.setCigarString("10M2N10M");
        Assert.assertFalse(CIGAR_IS_SUPPORTED.test(read), read.toString());
    }

    @DataProvider(name = "badCigars")
    public Object[][] badCigars() {
        return new Object[][]{
                {"2D4M"},               // starting with multiple deletions
                {"4M2D"},               // ending with multiple deletions
                {"3M1I1D"},             // adjacent indels AND ends in deletion
                {"1M1I1D2M"},           // adjacent indels I->D
                {"1M1D2I1M"},           // adjacent indels D->I
                {"1M1I2M1D"},           // ends in single deletion with insertion in the middle
                {"4M1D"},               // ends in single deletion
                {"1D4M"},               // starts with single deletion
                {"2M1D1D2M"},           // adjacent D's
                {"1M1I1I1M"},           // adjacent I's
                {"1H1D4M"},             // starting with deletion after H
                {"1S1D3M"},             // starting with deletion after S
                {"1H1S1D3M"},           // starting with deletion after HS
                {"4M1D1H"},             // ending with deletion before H
                {"3M1D1S"},             // ending with deletion before S
                {"3M1D1S1H"},           // ending with deletion before HS
                {"1H1S1H1M"},           // H in the middle, after S
                {"1M1H1S1M"},           // S in the middle, after H
                {"10M2H10M"},           // H in the middle
                {"10M2S10M"},           // S in the middle
                {"1S1H"},               // only clipping
                {"1S1S"},               // only clipping
                {"1H1S"},               // only clipping
                {"1H1H"},               // only clipping
                {"1S1H10M"},            // H in the middle
                {"1H1M1S1M"},           // H in the middle
                {"1M1H1S"},             // H in the middle
                {"1H1S10M2S10M1S1H"},    // deceiving S in the middle: HSMSMSH
                {"1H1S10M2H10M1S1H"},    // deceiving H in the middle
                {"1H1H2M"},                   //  (invalid according to htsjdk)
                {"1S20S10M"},                 //  (invalid according to htsjdk)
                {"1S1S1S1M"},                 //  (invalid according to htsjdk)
                {"1H1S10M10S1S30H"},          //  (invalid according to htsjdk)
                {"1H1S10M10S1S30H"},          //  (invalid according to htsjdk)
                {"1H20H10M"},                 //  (invalid according to htsjdk)
                {"1H1H10M10H30H"},            //  (invalid according to htsjdk)
                {"1H1H10M10H1H30H"},          //  (invalid according to htsjdk)
                {"1M1H2H"},                   //  (invalid according to htsjdk)
        };
    }
    @Test(dataProvider = "badCigars")
    public void testWonkyCigars (String cigarString) {
        SAMRecord read = ReadClipperTestUtils.makeReadFromCigar(cigarString);
        Assert.assertFalse(GOOD_CIGAR.test(read), read.getCigarString());
    }

    @Test
    public void testReadCigarLengthMismatch() {
        SAMRecord read = ReadClipperTestUtils.makeReadFromCigar("4M", 1);
        Assert.assertFalse(READLENGTH_EQUALS_CIGARLENGTH.test(read), read.getCigarString());
    }

    @Test
    public void testEmptyCigar(){
        SAMRecord read = ReadClipperTestUtils.makeReadFromCigar("");
        Assert.assertTrue(GOOD_CIGAR.test(read), read.getCigarString());
    }

    @DataProvider(name = "goodCigars")
    public Object[][] goodCigars() {
        return new Object[][]{
                {"1H1S10M10S30H"},
                {"1I9H"},
                {"1I1S8H"},
                {"1S1I1S7H"}
        };
    }

    @Test(dataProvider = "goodCigars")
    public void testGoodCigars (String cigarString) {
        SAMRecord read = ReadClipperTestUtils.makeReadFromCigar(cigarString);
        Assert.assertTrue(GOOD_CIGAR.test(read), read.getCigarString());
    }
    @Test
    public void testGoodCigarsUpToSize() {
        //Note: not using data providers here because it's super slow to print (many minutes vs few seconds).
        List<Cigar> cigarList = ReadClipperTestUtils.generateCigarList(10);
        for (Cigar cigar : cigarList) {
            SAMRecord read = ReadClipperTestUtils.makeReadFromCigar(cigar);
            Assert.assertTrue(GOOD_CIGAR.test(read), read.getCigarString());
        }
    }
}
