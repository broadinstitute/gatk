package org.broadinstitute.hellbender.engine.filters;

import htsjdk.samtools.*;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;
import org.broadinstitute.hellbender.testutils.ReadClipperTestUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.lang.reflect.Method;
import java.util.*;

import static org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary.*;

/**
 * Tests for the read filter library.
 */
public final class ReadFilterLibraryUnitTest {
    private static final int CHR_COUNT = 2;
    private static final int CHR_START = 1;
    private static final int CHR_SIZE = 1000;
    private static final int GROUP_COUNT = 5;

    private SAMFileHeader createHeaderWithReadGroups() {
        return ArtificialReadUtils.createArtificialSamHeaderWithGroups(CHR_COUNT, CHR_START, CHR_SIZE, GROUP_COUNT);
    }

    /**
     * Creates a read record.
     *
     * @param header header for the new record
     * @param cigar the new record CIGAR.
     * @param group the new record group index that must be in the range \
     *              [0,{@link #GROUP_COUNT})
     * @param reference the reference sequence index (0-based)
     * @param start the start position of the read alignment in the reference
     *              (1-based)
     * @return never <code>null</code>
     */
    private GATKRead createRead( final SAMFileHeader header, final Cigar cigar, final int group, final int reference, final int start ) {
        final GATKRead record = ArtificialReadUtils.createArtificialRead(header, cigar);
        record.setPosition(header.getSequence(reference).getSequenceName(), start);
        record.setReadGroup(header.getReadGroups().get(group).getReadGroupId());
        return record;
    }

    private GATKRead createRead( final SAMFileHeader header, final String cigarString ) {
        final Cigar cigar = TextCigarCodec.decode(cigarString);
        return createRead(header, cigar, 1, 0, 10);
    }

    private GATKRead simpleGoodRead( final SAMFileHeader header ) {
        final String cigarString = "101M";
        final Cigar cigar = TextCigarCodec.decode(cigarString);
        GATKRead read = createRead(header, cigar, 1, 0, 10);
        read.setMappingQuality(50);
        return read;
    }

    @Test
    public void testCheckSeqStored() {
        final GATKRead goodRead = ArtificialReadUtils.createArtificialRead(new byte[]{(byte) 'A'}, new byte[]{(byte) 'A'}, "1M");
        final GATKRead badRead = ArtificialReadUtils.createArtificialRead(new byte[]{}, new byte[]{}, "1M");
        badRead.setBases(new byte[0]);

        Assert.assertTrue(SEQ_IS_STORED.test(goodRead));
        Assert.assertFalse(SEQ_IS_STORED.test(badRead));
    }

    @Test(dataProvider = "UnsupportedCigarOperatorDataProvider")
    public void testCigarNOperatorFilter(String cigarString) {
        final SAMFileHeader header = createHeaderWithReadGroups();
        final ReadFilter filter = new WellformedReadFilter(header);
        final GATKRead read = createRead(header, cigarString);
        final boolean containsN = cigarString.contains("N");
        Assert.assertEquals(containsN, !filter.test(read), cigarString);
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

    @Test
    public void passesAllFilters() {
        final SAMFileHeader header = createHeaderWithReadGroups();
        final GATKRead read = simpleGoodRead(header);

        Assert.assertTrue(MAPPED.test(read), "MAPPED " + read.toString());
        Assert.assertTrue(NOT_SECONDARY_ALIGNMENT.test(read), "NOT_SECONDARY_ALIGNMENT " + read.toString());
        Assert.assertTrue(NOT_DUPLICATE.test(read), "NOT_DUPLICATE " + read.toString());
        Assert.assertTrue(PASSES_VENDOR_QUALITY_CHECK.test(read), "PASSES_VENDOR_QUALITY_CHECK " + read.toString());
        Assert.assertTrue(MAPPING_QUALITY_AVAILABLE.test(read), "MAPPING_QUALITY_AVAILABLE " + read.toString());
        Assert.assertTrue(MAPPING_QUALITY_NOT_ZERO.test(read), "MAPPING_QUALITY_NOT_ZERO " + read.toString());
        Assert.assertTrue(VALID_ALIGNMENT_START.test(read), "VALID_ALIGNMENT_START " + read.toString());
        Assert.assertTrue(VALID_ALIGNMENT_END.test(read), "VALID_ALIGNMENT_END " + read.toString());
        Assert.assertTrue(HAS_READ_GROUP.test(read), "HAS_READ_GROUP " + read.toString());
        Assert.assertTrue(HAS_MATCHING_BASES_AND_QUALS.test(read), "HAS_MATCHING_BASES_AND_QUALS " + read.toString());
        Assert.assertTrue(SEQ_IS_STORED.test(read), "SEQ_IS_STORED " + read.toString());
        Assert.assertTrue(CIGAR_CONTAINS_NO_N_OPERATOR.test(read), "CIGAR_CONTAINS_NO_N_OPERATOR " + read.toString());

        final WellformedReadFilter wellformed = new WellformedReadFilter(header);
        Assert.assertTrue(wellformed.test(read), "WELLFORMED " + read.toString());
    }

    @Test
    public void failsMAPPED_flag() {
        final SAMFileHeader header = createHeaderWithReadGroups();
        final GATKRead read = simpleGoodRead(header);
        read.setIsUnmapped();
        Assert.assertFalse(MAPPED.test(read), read.toString());
    }

    @Test
    public void failsNOT_DUPLICATE() {
        final SAMFileHeader header = createHeaderWithReadGroups();
        final GATKRead read = simpleGoodRead(header);
        read.setIsDuplicate(true);
        Assert.assertFalse(NOT_DUPLICATE.test(read), read.toString());
    }

    @Test
    public void failsPASSES_VENDOR_QUALITY_CHECK() {
        final SAMFileHeader header = createHeaderWithReadGroups();
        final GATKRead read = simpleGoodRead(header);
        read.setFailsVendorQualityCheck(true);
        Assert.assertFalse(PASSES_VENDOR_QUALITY_CHECK.test(read), read.toString());
    }

    @Test
    public void failsMAPPING_QUALITY_AVAILABLE() {
        final SAMFileHeader header = createHeaderWithReadGroups();
        final GATKRead read = simpleGoodRead(header);
        read.setMappingQuality(QualityUtils.MAPPING_QUALITY_UNAVAILABLE);
        Assert.assertFalse(MAPPING_QUALITY_AVAILABLE.test(read), read.toString());
    }

    @Test
    public void failsMAPPING_QUALITY_NOT_ZERO() {
        final SAMFileHeader header = createHeaderWithReadGroups();
        final GATKRead read = simpleGoodRead(header);
        read.setMappingQuality(0);
        Assert.assertFalse(MAPPING_QUALITY_NOT_ZERO.test(read), read.toString());
    }

    @Test
    public void VALID_ALIGNMENT_START_allows_unmapped_case1() {
        final SAMFileHeader header = createHeaderWithReadGroups();
        final GATKRead read = simpleGoodRead(header);
        read.setIsUnmapped();
        Assert.assertTrue(VALID_ALIGNMENT_START.test(read), "VALID_ALIGNMENT_START failed on an unmapped read");
    }

    @Test
    public void VALID_ALIGNMENT_START_allows_unmapped_case2() {
        final SAMFileHeader header = createHeaderWithReadGroups();
        final GATKRead read = simpleGoodRead(header);
        final SAMRecord samWithUnmappedPosition = read.convertToSAMRecord(header); // GATKRead interface will not let us set invalid alignment starts,
                                                                                   // so we need to convert to SAMRecord here
        samWithUnmappedPosition.setAlignmentStart(0);

        // Reads with a start position of 0 are considered unmapped in our Read interface
        GATKRead readWithUnmappedPosition = new SAMRecordToGATKReadAdapter(samWithUnmappedPosition);
        Assert.assertTrue(VALID_ALIGNMENT_START.test(readWithUnmappedPosition), "VALID_ALIGNMENT_START failed on an unmapped read (with start position == 0)");
    }

    @Test
    public void failsVALID_ALIGNMENT_START() {
        final SAMFileHeader header = createHeaderWithReadGroups();
        final GATKRead read = simpleGoodRead(header);
        final SAMRecord corruptSAM = read.convertToSAMRecord(header); // GATKRead interface will not let us set invalid alignment starts,
                                                                      // so we need to convert to SAMRecord here
        corruptSAM.setAlignmentStart(-1);
        Assert.assertFalse(VALID_ALIGNMENT_START.test(new SAMRecordToGATKReadAdapter(corruptSAM)), read.toString());
    }

    @Test
    public void passesValidAlignmentEndReadConsumesZeroRefBases() {
        final SAMFileHeader header = createHeaderWithReadGroups();
        final GATKRead read = simpleGoodRead(header);

        // Read with only soft-clips and insertions, consuming no reference bases
        read.setCigar("60S30I11S");
        Assert.assertEquals(read.getCigar().getReferenceLength(), 0, "read should consume no reference bases");

        Assert.assertTrue(VALID_ALIGNMENT_END.test(read), "Read consuming no reference bases should have passed VALID_ALIGNMENT_END");
    }

    @Test
    public void failsHAS_READ_GROUP() {
        final SAMFileHeader header = createHeaderWithReadGroups();
        final GATKRead read = simpleGoodRead(header);
        read.setReadGroup(null);
        Assert.assertFalse(HAS_READ_GROUP.test(read), read.toString());
    }

    @Test
    public void failsHAS_MATCHING_BASES_AND_QUALS() {
        final SAMFileHeader header = createHeaderWithReadGroups();
        final GATKRead read = simpleGoodRead(header);
        read.setBaseQualities(new byte[]{1,2,3});
        Assert.assertFalse(HAS_MATCHING_BASES_AND_QUALS.test(read), read.toString());
    }

    @Test
    public void failsSEQ_IS_STORED() {
        final SAMFileHeader header = createHeaderWithReadGroups();
        final GATKRead read = simpleGoodRead(header);
        read.setBases(new byte[0]);
        Assert.assertFalse(SEQ_IS_STORED.test(read), read.toString());
    }

    @Test
    public void failsCIGAR_CONTAINS_NO_N_OPERATOR() {
        final SAMFileHeader header = createHeaderWithReadGroups();
        final GATKRead read = simpleGoodRead(header);
        read.setCigar("10M2N10M");
        Assert.assertFalse(CIGAR_CONTAINS_NO_N_OPERATOR.test(read), read.toString());
    }

    @Test(dataProvider = "nonZeroReferenceLengthAlignmentFilterData")
    public void testNonZeroReferenceLengthAlignmentFilter(final String cigarString, final boolean expected) {
        final SAMFileHeader header = createHeaderWithReadGroups();
        final ReadFilter filter = ReadFilterLibrary.NON_ZERO_REFERENCE_LENGTH_ALIGNMENT;
        final GATKRead read = createRead(header, cigarString);
        Assert.assertEquals(filter.test(read), expected, cigarString);
    }

    @DataProvider (name = "nonZeroReferenceLengthAlignmentFilterData")
    public Object[][] nonZeroReferenceLengthAlignmentFilterData() {
        return new Object[][] {
                {"101M", true},
                {"1M", true},
                {"0M", false},
                {"1D", true},
                {"0D", false},
                {"0D0M", false},
                {"50I20S", false},
                {"50I1M20S", true},
                {"10I0M20S2H", false},
                {"10H10S50I50P10S10H",false},
                {"10D",true},
                {"1M2D1M10I10S", true}
        };
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
        GATKRead read = ReadClipperTestUtils.makeReadFromCigar(cigarString);
        Assert.assertFalse(GOOD_CIGAR.test(read), read.getCigar().toString());
    }

    @Test
    public void testReadCigarLengthMismatch() {
        GATKRead read = ReadClipperTestUtils.makeReadFromCigar("4M", 1);
        Assert.assertFalse(READLENGTH_EQUALS_CIGARLENGTH.test(read), read.getCigar().toString());
    }

    @Test
    public void testReadCigarLengthMismatch_wellFormed() {
        final SAMFileHeader header = createHeaderWithReadGroups();
        final GATKRead read = simpleGoodRead(header);

        //The read passes the filter here
        Assert.assertTrue(new WellformedReadFilter(header).test(read), read.getCigar().toString());

        //now we mess up the read by adding some bases and quals
        final int len = read.getLength();
        read.setBases(Utils.dupBytes((byte)'A', len + 1));
        read.setBaseQualities(Utils.dupBytes((byte)60, len + 1));

        //And now the read fails the filter
        Assert.assertNotEquals(read.getCigar().getReadLength(), read.getLength());
        Assert.assertFalse(new WellformedReadFilter(header).test(read), read.getCigar().toString());
    }

    @Test
    public void testEmptyCigar(){
        GATKRead read = ReadClipperTestUtils.makeReadFromCigar("");
        Assert.assertTrue(GOOD_CIGAR.test(read), read.getCigar().toString());
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
        GATKRead read = ReadClipperTestUtils.makeReadFromCigar(cigarString);
        Assert.assertTrue(GOOD_CIGAR.test(read), read.getCigar().toString());
    }
    @Test
    public void testGoodCigarsUpToSize() {
        //Note: not using data providers here because it's super slow to print (many minutes vs few seconds).
        List<Cigar> cigarList = ReadClipperTestUtils.generateCigarList(10);
        for (Cigar cigar : cigarList) {
            GATKRead read = ReadClipperTestUtils.makeReadFromCigar(cigar);
            Assert.assertTrue(GOOD_CIGAR.test(read), read.getCigar().toString());
        }
    }

    @Test
    public void testAlignmentAgreesWithHeader() {
        final SAMFileHeader header = createHeaderWithReadGroups();
        final AlignmentAgreesWithHeaderReadFilter filter = new AlignmentAgreesWithHeaderReadFilter(header);
        final GATKRead read = simpleGoodRead(header);

        read.setPosition("BAD_CONTIG", 1);
        Assert.assertFalse(filter.test(read), "AlignmentAgreesWithHeader read filter should have failed on read with bad contig");

        read.setPosition(Integer.toString(CHR_START), CHR_SIZE + 1);
        Assert.assertFalse(filter.test(read), "AlignmentAgreesWithHeader read filter should have failed on read with start position past the end of its contig");
    }

    @Test
    public void testLibraryReadFilter() {
        final SAMFileHeader header = createHeaderWithReadGroups();
        final GATKRead read = simpleGoodRead(header);
        final LibraryReadFilter f = new LibraryReadFilter(header);

        final String foo = "Foo";
        header.getReadGroup(read.getReadGroup()).setLibrary(foo);

        Assert.assertFalse(f.test(read), read.toString());//fail
        f.libraryToKeep = Collections.singleton(foo);
        Assert.assertTrue(f.test(read), read.toString());//pass
        f.libraryToKeep = new HashSet<>(Arrays.asList("A", "B"));
        Assert.assertFalse(f.test(read), read.toString());
        f.libraryToKeep.add(foo);
        Assert.assertTrue(f.test(read), read.toString());
    }

    @Test
    public void testMappingQualityFilter() {
        final SAMFileHeader header = createHeaderWithReadGroups();
        final GATKRead read = simpleGoodRead(header);

        MappingQualityReadFilter f = new MappingQualityReadFilter(17);
        read.setMappingQuality(11);
        Assert.assertFalse(f.test(read), read.toString());//fail

        f = new MappingQualityReadFilter(9);
        Assert.assertTrue(f.test(read), read.toString());//pass

        // with maximum mapping quality
        f = new MappingQualityReadFilter(1, 10);
        Assert.assertFalse(f.test(read), read.toString());//fail

        f = new MappingQualityReadFilter(9, 12);
        Assert.assertTrue(f.test(read), read.toString());//pass

        // limit range to the same mapping quality
        f = new MappingQualityReadFilter(11, 11);
        Assert.assertTrue(f.test(read), read.toString());//pass

        // limit range to lower/higher mapping quality
        f = new MappingQualityReadFilter(10, 10);
        Assert.assertFalse(f.test(read), read.toString());//fail

        f = new MappingQualityReadFilter(12, 12);
        Assert.assertFalse(f.test(read), read.toString());//fail
    }

    @Test
    public void testMaxInsertSizeFilter() {
        final SAMFileHeader header = createHeaderWithReadGroups();
        final GATKRead pairedRead = simpleGoodRead(header);
        final GATKRead unpairedRead = simpleGoodRead(header);
        pairedRead.setIsPaired(true);

        final FragmentLengthReadFilter f = new FragmentLengthReadFilter();

        pairedRead.setFragmentLength(150);
        unpairedRead.setFragmentLength(150);

        f.maxFragmentLength = 180;
        Assert.assertTrue(f.test(pairedRead), pairedRead.toString());//pass
        Assert.assertTrue(f.test(unpairedRead), pairedRead.toString());//pass

        f.maxFragmentLength = 90;
        Assert.assertFalse(f.test(pairedRead), pairedRead.toString());//fail
        Assert.assertTrue(f.test(unpairedRead), pairedRead.toString());//pass

        pairedRead.setFragmentLength(-150);

        f.maxFragmentLength = 180;
        Assert.assertTrue(f.test(pairedRead), pairedRead.toString());//pass
        Assert.assertTrue(f.test(unpairedRead), pairedRead.toString());//pass

        f.maxFragmentLength = 90;
        Assert.assertFalse(f.test(pairedRead), pairedRead.toString());//fail
        Assert.assertTrue(f.test(unpairedRead), pairedRead.toString());//pass
    }

    @Test
    public void testPlatformFilter() {
        final SAMFileHeader header = createHeaderWithReadGroups();
        final GATKRead read = simpleGoodRead(header);
        PlatformReadFilter f = new PlatformReadFilter(header);

        f.PLFilterNames = new LinkedHashSet<>(Arrays.asList("PL1", "PL2"));
        header.getReadGroup(read.getReadGroup()).setPlatform("PL1");
        Assert.assertTrue(f.test(read), read.toString());//pass

        header.getReadGroup(read.getReadGroup()).setPlatform(null);
        Assert.assertFalse(f.test(read), read.toString());//fail - no match

        header.getReadGroup(read.getReadGroup()).setPlatform("prefix pl1 suffix");  //not exact matching
        Assert.assertTrue(f.test(read), read.toString());//pass

        f.PLFilterNames = new LinkedHashSet<>(Arrays.asList("Fred"));
        header.getReadGroup(read.getReadGroup()).setPlatform("PL1");
        Assert.assertFalse(f.test(read), read.toString());//fail
    }

    @Test
    public void testReadLengthFilter() {
        final SAMFileHeader header = createHeaderWithReadGroups();
        final GATKRead read = simpleGoodRead(header);
        ReadLengthReadFilter f = new ReadLengthReadFilter();
        f.minReadLength = 10;
        f.maxReadLength = 20;

        read.setBases(new byte[5]);
        Assert.assertFalse(f.test(read), read.toString());//fail

        read.setBases(new byte[10]);
        Assert.assertTrue(f.test(read), read.toString());//pass

        read.setBases(new byte[15]);
        Assert.assertTrue(f.test(read), read.toString());//pass

        read.setBases(new byte[20]);
        Assert.assertTrue(f.test(read), read.toString());//pass

        read.setBases(new byte[25]);
        Assert.assertFalse(f.test(read), read.toString());//fail
    }

    @Test
    public void testReadNameFilter() {
        final SAMFileHeader header = createHeaderWithReadGroups();
        final GATKRead read = simpleGoodRead(header);
        ReadNameReadFilter f = new ReadNameReadFilter();

        final String fred= "fred";
        f.readName = fred;
        read.setName(fred);
        Assert.assertTrue(f.test(read), read.toString());//pass

        read.setName(fred.toUpperCase());
        Assert.assertFalse(f.test(read), read.toString());//fail
    }

    @Test
    public void testReadStrandFilter() {
        final SAMFileHeader header = createHeaderWithReadGroups();
        final GATKRead read = simpleGoodRead(header);
        ReadStrandFilter f = new ReadStrandFilter();

        f.keepOnlyReverse = false;
        read.setIsReverseStrand(false);
        Assert.assertTrue(f.test(read), read.toString());//pass

        read.setIsReverseStrand(true);
        Assert.assertFalse(f.test(read), read.toString());//fail

        f.keepOnlyReverse = true;
        read.setIsReverseStrand(false);
        Assert.assertFalse(f.test(read), read.toString());//fail

        read.setIsReverseStrand(true);
        Assert.assertTrue(f.test(read), read.toString());//pass
    }

    @Test
    public void testSampleFilter() {
        final SAMFileHeader header = createHeaderWithReadGroups();
        final GATKRead read = simpleGoodRead(header);
        SampleReadFilter f = new SampleReadFilter(header);

        final String fred = "fred";
        f.samplesToKeep = Collections.emptySet();
        header.getReadGroup(read.getReadGroup()).setSample(fred);
        Assert.assertFalse(f.test(read), read.toString());//fail

        f.samplesToKeep = Collections.singleton(fred);
        header.getReadGroup(read.getReadGroup()).setSample(fred);
        Assert.assertTrue(f.test(read), read.toString());//pass

        f.samplesToKeep = Collections.singleton(fred);
        header.getReadGroup(read.getReadGroup()).setSample(fred + "suffix");
        Assert.assertFalse(f.test(read), read.toString());//fail - exact matching

        f.samplesToKeep = Collections.singleton(fred);
        header.getReadGroup(read.getReadGroup()).setSample(fred.toUpperCase());
        Assert.assertFalse(f.test(read), read.toString());//fail - case sensitive matching

        f.samplesToKeep = new LinkedHashSet<>(Arrays.asList(fred, "bozo"));
        header.getReadGroup(read.getReadGroup()).setSample(fred);
        Assert.assertTrue(f.test(read), read.toString());//pass
    }

    @Test
    public void testSingleReadGroupFilter() {
        final SAMFileHeader header = createHeaderWithReadGroups();
        final GATKRead read = simpleGoodRead(header);
        final String fred = "fred";

        ReadGroupReadFilter f = new ReadGroupReadFilter();

        f.readGroup = "";
        read.setReadGroup(fred);
        Assert.assertFalse(f.test(read), read.toString());//fail

        f.readGroup = fred;
        Assert.assertTrue(f.test(read), read.toString());//pass

        f.readGroup = fred;
        read.setReadGroup(fred + "suffix");
        Assert.assertFalse(f.test(read), read.toString());//fail - exact matching

        f.readGroup = fred;
        read.setReadGroup(fred.toUpperCase());
        Assert.assertFalse(f.test(read), read.toString());//fail - case sensitive matching
    }

    @Test
    public void testPlatformUnitFilter() {
        final SAMFileHeader header = createHeaderWithReadGroups();
        final PlatformUnitReadFilter f = new PlatformUnitReadFilter(header);
        final GATKRead read = simpleGoodRead(header);
        final String fred = "fred";
        header.getReadGroup(read.getReadGroup()).setPlatformUnit(fred);

        f.blackListedLanes = Collections.emptySet();
        header.getReadGroup(read.getReadGroup()).setPlatformUnit(fred);
        Assert.assertTrue(f.test(read), read.toString());//pass - no blacklist

        f.blackListedLanes = Collections.singleton(fred);
        header.getReadGroup(read.getReadGroup()).setPlatformUnit(fred);
        Assert.assertFalse(f.test(read), read.toString());//fail

        f.blackListedLanes = Collections.singleton(fred);
        header.getReadGroup(read.getReadGroup()).setPlatformUnit(fred + "suffix");
        Assert.assertTrue(f.test(read), read.toString());//pass - exact matching

        f.blackListedLanes = Collections.singleton(fred);
        header.getReadGroup(read.getReadGroup()).setPlatformUnit(fred.toUpperCase());
        Assert.assertTrue(f.test(read), read.toString());//pass - case sensitive matching

        f.blackListedLanes = new LinkedHashSet<>(Arrays.asList(fred, "bozo"));
        header.getReadGroup(read.getReadGroup()).setPlatformUnit(fred);
        Assert.assertFalse(f.test(read), read.toString());//fail

        f.blackListedLanes = Collections.singleton(fred);
        read.setAttribute(SAMTag.PU.name(), fred);
        Assert.assertFalse(f.test(read), read.toString());//fail - match
    }

    @DataProvider(name = "MateOnSameContigOrNoMappedMateTestData")
    public Object[][] mateOnSameContigOrNoMappedMateTestData() {
        final SAMFileHeader header = createHeaderWithReadGroups();

        final GATKRead unpairedRead = simpleGoodRead(header);
        unpairedRead.setIsPaired(false);

        final GATKRead pairedReadWithUnmappedMate = simpleGoodRead(header);
        pairedReadWithUnmappedMate.setIsPaired(true);
        pairedReadWithUnmappedMate.setMateIsUnmapped();

        final GATKRead pairedReadWithMappedMateDifferentContig = simpleGoodRead(header);
        pairedReadWithMappedMateDifferentContig.setIsPaired(true);
        pairedReadWithMappedMateDifferentContig.setMatePosition("2", 1);
        pairedReadWithMappedMateDifferentContig.setPosition("1", 1);

        final GATKRead pairedReadWithMappedMateSameContig = simpleGoodRead(header);
        pairedReadWithMappedMateSameContig.setIsPaired(true);
        pairedReadWithMappedMateSameContig.setMatePosition("1", 100);
        pairedReadWithMappedMateSameContig.setPosition("1", 1);

        final GATKRead unmappedReadWithMappedMate = simpleGoodRead(header);
        unmappedReadWithMappedMate.setIsUnmapped();
        unmappedReadWithMappedMate.setIsPaired(true);
        unmappedReadWithMappedMate.setMatePosition("1", 100);

        final GATKRead unmappedReadWithUnmappedMate = simpleGoodRead(header);
        unmappedReadWithUnmappedMate.setIsUnmapped();
        unmappedReadWithUnmappedMate.setIsPaired(true);
        unmappedReadWithUnmappedMate.setMateIsUnmapped();

        final GATKRead unmappedUnpairedRead = simpleGoodRead(header);
        unmappedUnpairedRead.setIsUnmapped();
        unmappedUnpairedRead.setIsPaired(false);

        return new Object[][] {
                { unpairedRead, true },
                { pairedReadWithUnmappedMate, true },
                { pairedReadWithMappedMateDifferentContig, false },
                { pairedReadWithMappedMateSameContig, true },
                { unmappedReadWithMappedMate, false },
                { unmappedReadWithUnmappedMate, true },
                { unmappedUnpairedRead, true }
        };
    }

    @Test(dataProvider = "MateOnSameContigOrNoMappedMateTestData")
    public void testMateOnSameContigOrMateNotMapped( final GATKRead read, final boolean expectedFilterResult ) {
        final boolean actualFilterResult = MATE_ON_SAME_CONTIG_OR_NO_MAPPED_MATE.test(read);
        Assert.assertEquals(actualFilterResult, expectedFilterResult);
    }

    @Test
    public void testFirstOfPair() {
        final SAMFileHeader header = createHeaderWithReadGroups();
        final GATKRead read = simpleGoodRead(header);

        read.setIsPaired(false);
        Assert.assertFalse(FIRST_OF_PAIR.test(read), "FIRST_OF_PAIR " + read.toString());
        read.setIsPaired(true);
        Assert.assertFalse(FIRST_OF_PAIR.test(read), "FIRST_OF_PAIR " + read.toString());
        read.setIsFirstOfPair();
        Assert.assertTrue(FIRST_OF_PAIR.test(read), "FIRST_OF_PAIR " + read.toString());
    }

    @Test
    public void testSecondOfPair() {
        final SAMFileHeader header = createHeaderWithReadGroups();
        final GATKRead read = simpleGoodRead(header);

        Assert.assertFalse(SECOND_OF_PAIR.test(read), "SECOND_OF_PAIR " + read.toString());
        read.setIsPaired(true);
        read.setIsSecondOfPair();
        Assert.assertTrue(SECOND_OF_PAIR.test(read), "SECOND_OF_PAIR " + read.toString());
    }

    @Test
    public void testMateDifferentStrandReadFilter() {
        final SAMFileHeader header = createHeaderWithReadGroups();
        final GATKRead pairedRead = simpleGoodRead(header);
        pairedRead.setIsPaired(true);

        Assert.assertFalse(MATE_DIFFERENT_STRAND.test(pairedRead), "MATE_DIFFERENT_STRAND " + pairedRead.toString());
        pairedRead.setIsReverseStrand(false);
        Assert.assertFalse(MATE_DIFFERENT_STRAND.test(pairedRead), "MATE_DIFFERENT_STRAND " + pairedRead.toString());
        pairedRead.setMatePosition("0", 1);
        Assert.assertFalse(MATE_DIFFERENT_STRAND.test(pairedRead), "MATE_DIFFERENT_STRAND " + pairedRead.toString());

        pairedRead.setMateIsReverseStrand(true);
        Assert.assertTrue(MATE_DIFFERENT_STRAND.test(pairedRead), "MATE_DIFFERENT_STRAND " + pairedRead.toString());
    }

    @Test
    public void testNonZeroFragmentLengthReadFilter() {
        final SAMFileHeader header = createHeaderWithReadGroups();
        final GATKRead read = simpleGoodRead(header);

        read.setFragmentLength(0);
        Assert.assertFalse(NONZERO_FRAGMENT_LENGTH_READ_FILTER.test(read), "NONZERO_FRAGMENT_LENGTH_READ_FILTER " + read.toString());
        read.setFragmentLength(9);
        Assert.assertTrue(NONZERO_FRAGMENT_LENGTH_READ_FILTER.test(read), "NONZERO_FRAGMENT_LENGTH_READ_FILTER " + read.toString());
    }

    @Test
    public void testNotSecondaryAlignmentReadFilter() {
        final SAMFileHeader header = createHeaderWithReadGroups();
        final GATKRead read = simpleGoodRead(header);

        read.setIsSecondaryAlignment(false);
        Assert.assertTrue(NOT_SECONDARY_ALIGNMENT.test(read), "NOT_SECONDARY_ALIGNMENT " + read.toString());
        read.setIsSecondaryAlignment(true);
        Assert.assertFalse(NOT_SECONDARY_ALIGNMENT.test(read), "NOT_SECONDARY_ALIGNMENT " + read.toString());
    }

    @Test
    public void testNotSupplementaryAlignmentReadFilter() {
        final SAMFileHeader header = createHeaderWithReadGroups();
        final GATKRead read = simpleGoodRead(header);

        read.setIsSupplementaryAlignment(false);
        Assert.assertTrue(NOT_SUPPLEMENTARY_ALIGNMENT.test(read), "NOT_SUPPLEMENTARY_ALIGNMENT " + read.toString());
        read.setIsSupplementaryAlignment(true);
        Assert.assertFalse(NOT_SUPPLEMENTARY_ALIGNMENT.test(read), "NOT_SUPPLEMENTARY_ALIGNMENT " + read.toString());
    }

    @Test
    public void testPairedReadFilter() {
        final SAMFileHeader header = createHeaderWithReadGroups();
        final GATKRead read = simpleGoodRead(header);
        read.setIsPaired(false);
        Assert.assertFalse(PAIRED.test(read), "PAIRED " + read.toString());
        read.setIsPaired(true);
        Assert.assertTrue(PAIRED.test(read), "PAIRED " + read.toString());
    }

    @Test
    public void testProperlyPairedReadFilter() {
        final SAMFileHeader header = createHeaderWithReadGroups();
        final GATKRead read = simpleGoodRead(header);
        Assert.assertFalse(PROPERLY_PAIRED.test(read), "PROPERLY_PAIRED " + read.toString());
        read.setIsPaired(true);
        read.setIsProperlyPaired(true);
        Assert.assertTrue(PROPERLY_PAIRED.test(read), "PROPERLY_PAIRED " + read.toString());
    }

    @Test
    public void testPrimaryAlignmentReadFilter() {
        // simple primary read (pass)
        final GATKRead read = simpleGoodRead(createHeaderWithReadGroups());
        Assert.assertTrue(ReadFilterLibrary.PRIMARY_LINE.test(read));

        // supplementary read (filter out)
        read.setIsSupplementaryAlignment(true);
        read.setIsSecondaryAlignment(false);
        Assert.assertFalse(ReadFilterLibrary.PRIMARY_LINE.test(read));

        // only secondary (filter out)
        read.setIsSupplementaryAlignment(false);
        read.setIsSecondaryAlignment(true);
        Assert.assertFalse(ReadFilterLibrary.PRIMARY_LINE.test(read));

        // both supplementary and secondary (filter out)
        read.setIsSupplementaryAlignment(true);
        read.setIsSecondaryAlignment(true);
        Assert.assertFalse(ReadFilterLibrary.PRIMARY_LINE.test(read));

    }

}
