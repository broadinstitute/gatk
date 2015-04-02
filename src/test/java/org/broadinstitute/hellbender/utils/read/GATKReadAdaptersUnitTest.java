package org.broadinstitute.hellbender.utils.read;


import com.google.api.services.genomics.model.Position;
import com.google.api.services.genomics.model.Read;
import htsjdk.samtools.*;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.List;

public class GATKReadAdaptersUnitTest extends BaseTest {

    public static final String BASIC_READ_NAME = "basic_read";
    public static final String BASIC_READ_CONTIG = "1";
    public static final int BASIC_READ_START = 5;
    public static final int BASIC_READ_MAPPING_QUALITY = 20;
    public static final byte[] BASIC_READ_BASES = {'A', 'C', 'G', 'T'};
    public static final byte[] BASIC_READ_BASE_QUALITIES = {30, 40, 30, 50};
    public static final String BASIC_READ_CIGAR = "1M1I2M";
    public static final String BASIC_READ_GROUP = "FOO";
    public static final int BASIC_READ_END = BASIC_READ_START + TextCigarCodec.decode(BASIC_READ_CIGAR).getReferenceLength() - 1;
    public static final String BASIC_READ_MATE_CONTIG = "1";
    public static final int BASIC_READ_MATE_START = 125;

    private static GATKRead basicReadBackedBySam() {
        return new SAMRecordToGATKReadAdapter(basicSAMRecord());
    }

    private static GATKRead basicReadBackedByGoogle() {
        return new GoogleGenomicsReadToGATKReadAdapter(basicGoogleGenomicsRead());
    }

    private static SAMRecord basicSAMRecord() {
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader(2, 1, 1000000);
        header.addReadGroup(new SAMReadGroupRecord("FOO"));

        final SAMRecord read = ArtificialReadUtils.createArtificialSAMRecord(
                header,
                BASIC_READ_NAME,
                header.getSequenceIndex(BASIC_READ_CONTIG),
                BASIC_READ_START,
                BASIC_READ_BASES,
                BASIC_READ_BASE_QUALITIES,
                BASIC_READ_CIGAR
        );
        read.setAttribute(SAMTag.RG.name(), BASIC_READ_GROUP);
        read.setMappingQuality(BASIC_READ_MAPPING_QUALITY);
        read.setMateReferenceName(BASIC_READ_MATE_CONTIG);
        read.setMateAlignmentStart(BASIC_READ_MATE_START);
        read.setReadPairedFlag(true);
        read.setFirstOfPairFlag(true);

        return read;
    }

    private static Read basicGoogleGenomicsRead() {
        final Read read = ArtificialReadUtils.createArtificialGoogleGenomicsRead(
                BASIC_READ_NAME,
                BASIC_READ_CONTIG,
                BASIC_READ_START,
                BASIC_READ_BASES,
                BASIC_READ_BASE_QUALITIES,
                BASIC_READ_CIGAR
        );
        read.setReadGroupId(BASIC_READ_GROUP);
        read.getAlignment().getPosition().setReverseStrand(false);
        read.getAlignment().setMappingQuality(BASIC_READ_MAPPING_QUALITY);
        read.setNextMatePosition(new Position());
        read.getNextMatePosition().setReferenceName(BASIC_READ_MATE_CONTIG);
        read.getNextMatePosition().setPosition((long) BASIC_READ_MATE_START - 1);
        read.getNextMatePosition().setReverseStrand(false);
        read.setNumberReads(2);
        read.setReadNumber(0);
        read.setProperPlacement(false);

        return read;
    }

    private static List<GATKRead> getUnmappedReads() {
        List<GATKRead> unmappedReads = new ArrayList<>();

        final SAMRecord unmappedFlagSam = basicSAMRecord();
        unmappedFlagSam.setReadUnmappedFlag(true);
        unmappedReads.add(new SAMRecordToGATKReadAdapter(unmappedFlagSam));

        final SAMRecord unmappedContigSam = basicSAMRecord();
        unmappedContigSam.setReferenceName(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME);
        unmappedReads.add(new SAMRecordToGATKReadAdapter(unmappedContigSam));

        final SAMRecord noAlignmentStartSam = basicSAMRecord();
        noAlignmentStartSam.setAlignmentStart(SAMRecord.NO_ALIGNMENT_START);
        unmappedReads.add(new SAMRecordToGATKReadAdapter(noAlignmentStartSam));

        final Read noAlignmentGoogleRead = basicGoogleGenomicsRead();
        noAlignmentGoogleRead.setAlignment(null);
        unmappedReads.add(new GoogleGenomicsReadToGATKReadAdapter(noAlignmentGoogleRead));

        final Read noPositionGoogleRead = basicGoogleGenomicsRead();
        noPositionGoogleRead.getAlignment().setPosition(null);
        unmappedReads.add(new GoogleGenomicsReadToGATKReadAdapter(noPositionGoogleRead));

        final Read noContigGoogleRead = basicGoogleGenomicsRead();
        noContigGoogleRead.getAlignment().getPosition().setReferenceName(null);
        unmappedReads.add(new GoogleGenomicsReadToGATKReadAdapter(noContigGoogleRead));

        final Read starContigGoogleRead = basicGoogleGenomicsRead();
        starContigGoogleRead.getAlignment().getPosition().setReferenceName(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME);
        unmappedReads.add(new GoogleGenomicsReadToGATKReadAdapter(starContigGoogleRead));

        final Read noStartGoogleRead = basicGoogleGenomicsRead();
        noStartGoogleRead.getAlignment().getPosition().setPosition(-1l);
        unmappedReads.add(new GoogleGenomicsReadToGATKReadAdapter(noStartGoogleRead));

        return unmappedReads;
    }

    @DataProvider(name = "GetAndSetPositionData")
    public Object[][] getAndSetPositionData() {
        List<Object[]> testCases = new ArrayList<>();

        testCases.add(new Object[]{basicReadBackedBySam(), BASIC_READ_CONTIG, BASIC_READ_START, BASIC_READ_END});
        testCases.add(new Object[]{basicReadBackedByGoogle(), BASIC_READ_CONTIG, BASIC_READ_START, BASIC_READ_END});

        for ( GATKRead unmappedRead : getUnmappedReads() ) {
            testCases.add(new Object[]{unmappedRead, null, ReadConstants.UNSET_POSITION, ReadConstants.UNSET_POSITION});
        }

        return testCases.toArray(new Object[][]{});
    }

    @Test(dataProvider = "GetAndSetPositionData")
    public void testGetAndSetPosition( final GATKRead read, final String expectedContig, final int expectedStart, final int expectedEnd ) {
        Assert.assertEquals(read.getContig(), expectedContig, "Wrong contig for read");
        Assert.assertEquals(read.getStart(), expectedStart, "Wrong start position for read");
        Assert.assertEquals(read.getEnd(), expectedEnd, "Wrong end position for read");

        read.setPosition("2", 17);
        Assert.assertFalse(read.isUnmapped(), "Read should be mapped after setPosition() call");
        Assert.assertEquals(read.getContig(), "2", "Wrong contig after setPosition()");
        Assert.assertEquals(read.getStart(), 17, "Wrong start position after setPosition()");
    }

    @DataProvider(name = "InvalidSetPositionData")
    public Object[][] invalidSetPositionData() {
        return new Object[][]{
                { null, 1 },
                { SAMRecord.NO_ALIGNMENT_REFERENCE_NAME, 1 },
                { "1", SAMRecord.NO_ALIGNMENT_START }
        };
    }

    @Test(dataProvider = "InvalidSetPositionData", expectedExceptions = IllegalArgumentException.class)
    public void testSetInvalidPositionOnSamBackedRead( final String contig, final int start ) {
        final GATKRead read = new SAMRecordToGATKReadAdapter(basicSAMRecord());
        read.setPosition(contig, start);
    }

    @Test(dataProvider = "InvalidSetPositionData", expectedExceptions = IllegalArgumentException.class)
    public void testSetInvalidPositionOnGoogleBackedRead( final String contig, final int start ) {
        final GATKRead read = new GoogleGenomicsReadToGATKReadAdapter(basicGoogleGenomicsRead());
        read.setPosition(contig, start);
    }

    @DataProvider(name = "GetAndSetNameData")
    public Object[][] getAndSetNameData() {
        final SAMRecord namelessSam = basicSAMRecord();
        namelessSam.setReadName(null);

        final Read namelessGoogleRead = basicGoogleGenomicsRead();
        namelessGoogleRead.setFragmentName(null);

        return new Object[][]{
                { basicReadBackedBySam(), BASIC_READ_NAME },
                { basicReadBackedByGoogle(), BASIC_READ_NAME },
                { new SAMRecordToGATKReadAdapter(namelessSam), null },
                { new GoogleGenomicsReadToGATKReadAdapter(namelessGoogleRead), null }
        };
    }

    @Test(dataProvider = "GetAndSetNameData")
    public void testGetAndSetName( final GATKRead read, final String expectedName ) {
        Assert.assertEquals(read.getName(), expectedName, "Wrong name for read");

        read.setName("NEWNAME");
        Assert.assertEquals(read.getName(), "NEWNAME", "Wrong name for read after setName()");

        read.setName(null);
        Assert.assertEquals(read.getName(), null, "Wrong name for read after setName(null)");
    }

    @DataProvider(name = "GetLengthData")
    public Object[][] getLengthData() {
        final SAMRecord baselessSam = basicSAMRecord();
        baselessSam.setReadBases(SAMRecord.NULL_SEQUENCE);

        final Read baselessGoogleRead = basicGoogleGenomicsRead();
        baselessGoogleRead.setAlignedSequence(null);

        return new Object[][]{
                { basicReadBackedBySam(), BASIC_READ_BASES.length },
                { basicReadBackedByGoogle(), BASIC_READ_BASES.length },
                { new SAMRecordToGATKReadAdapter(baselessSam), 0 },
                { new GoogleGenomicsReadToGATKReadAdapter(baselessGoogleRead), 0 }
        };
    }

    @Test(dataProvider = "GetLengthData")
    public void testGetLength( final GATKRead read, final int expectedLength ) {
        Assert.assertEquals(read.getLength(), expectedLength, "Wrong length for read");
    }

    @DataProvider(name = "GetUnclippedStartAndEndData")
    public Object[][] getUnclippedStartAndEndData() {
        final SAMRecord softClippedSam = basicSAMRecord();
        softClippedSam.setCigarString("1S2M1S");

        final SAMRecord hardClippedSam = basicSAMRecord();
        hardClippedSam.setCigarString("3H2M2H");

        final Read softClippedGoogleRead = basicGoogleGenomicsRead();
        softClippedGoogleRead.getAlignment().setCigar(CigarConversionUtils.convertSAMCigarToCigarUnitList(TextCigarCodec.decode("1S2M1S")));

        final Read hardClippedGoogleRead = basicGoogleGenomicsRead();
        hardClippedGoogleRead.getAlignment().setCigar(CigarConversionUtils.convertSAMCigarToCigarUnitList(TextCigarCodec.decode("3H2M2H")));

        return new Object[][]{
                { new SAMRecordToGATKReadAdapter(softClippedSam), BASIC_READ_START - 1, BASIC_READ_START + 2 },
                { new SAMRecordToGATKReadAdapter(hardClippedSam), BASIC_READ_START - 3, BASIC_READ_START + 3 },
                { new GoogleGenomicsReadToGATKReadAdapter(softClippedGoogleRead), BASIC_READ_START - 1, BASIC_READ_START + 2 },
                { new GoogleGenomicsReadToGATKReadAdapter(hardClippedGoogleRead), BASIC_READ_START - 3, BASIC_READ_START + 3 }
        };
    }

    @Test(dataProvider = "GetUnclippedStartAndEndData")
    public void testGetAndSetUnclippedStartAndEnd( final GATKRead read, final int expectedUnclippedStart, final int expectedUnclippedEnd ) {
        Assert.assertEquals(read.getUnclippedStart(), expectedUnclippedStart, "Wrong unclipped start for read");
        Assert.assertEquals(read.getUnclippedEnd(), expectedUnclippedEnd, "Wrong unclipped end for read");
    }

    @DataProvider(name = "GetAndSetMatePositionData")
    public Object[][] getAndSetMatePositionData() {
        final SAMRecord samWithUnmappedMate = basicSAMRecord();
        samWithUnmappedMate.setMateUnmappedFlag(true);

        final Read googleReadWithUnmappedMate = basicGoogleGenomicsRead();
        googleReadWithUnmappedMate.setNextMatePosition(null);

        return new Object[][]{
                { basicReadBackedBySam(), BASIC_READ_MATE_CONTIG, BASIC_READ_MATE_START },
                { basicReadBackedByGoogle(), BASIC_READ_MATE_CONTIG, BASIC_READ_MATE_START },
                { new SAMRecordToGATKReadAdapter(samWithUnmappedMate), null, ReadConstants.UNSET_POSITION },
                { new GoogleGenomicsReadToGATKReadAdapter(googleReadWithUnmappedMate), null, ReadConstants.UNSET_POSITION }
        };
    }

    @Test(dataProvider = "GetAndSetMatePositionData")
    public void testGetAndSetMatePosition( final GATKRead read, final String expectedMateContig, final int expectedMateStart ) {
        Assert.assertEquals(read.getMateContig(), expectedMateContig, "Wrong contig for mate");
        Assert.assertEquals(read.getMateStart(), expectedMateStart, "Wrong start position for mate");

        read.setMatePosition("2", 52);
        Assert.assertFalse(read.mateIsUnmapped(), "Mate should be mapped after setMatePosition()");
        Assert.assertEquals(read.getMateContig(), "2", "Wrong contig for mate after setMatePosition()");
        Assert.assertEquals(read.getMateStart(), 52, "Wrong start position for mate after setMatePosition()");

        read.setMatePosition(new SimpleInterval("1", 40, 40));
        Assert.assertEquals(read.getMateContig(), "1", "Wrong contig for mate after setMatePosition()");
        Assert.assertEquals(read.getMateStart(), 40, "Wrong start position for mate after setMatePosition()");

        // Setting mate position should have the additional effect of marking the read as paired
        read.setIsPaired(false);
        read.setMatePosition("1", 1);
        Assert.assertTrue(read.isPaired(), "Read should be marked as paired after setMatePosition()");

        read.setIsPaired(false);
        read.setMatePosition(new SimpleInterval("2", 1, 1));
        Assert.assertTrue(read.isPaired(), "Read should be marked as paired after setMatePosition()");
    }

    @Test(dataProvider = "InvalidSetPositionData", expectedExceptions = IllegalArgumentException.class)
    public void testSetInvalidMatePositionOnSamBackedRead( final String contig, final int start ) {
        final GATKRead read = new SAMRecordToGATKReadAdapter(basicSAMRecord());
        read.setMatePosition(contig, start);
    }

    @Test(dataProvider = "InvalidSetPositionData", expectedExceptions = IllegalArgumentException.class)
    public void testSetInvalidMatePositionOnGoogleBackedRead( final String contig, final int start ) {
        final GATKRead read = new GoogleGenomicsReadToGATKReadAdapter(basicGoogleGenomicsRead());
        read.setMatePosition(contig, start);
    }

    @DataProvider(name = "GetAndSetFragmentLengthData")
    public Object[][] getAndSetFragmentLengthData() {
        final SAMRecord samWithISize = basicSAMRecord();
        samWithISize.setInferredInsertSize(120);

        final Read googleReadWithFragmentLength = basicGoogleGenomicsRead();
        googleReadWithFragmentLength.setFragmentLength(120);

        return new Object[][]{
                { basicReadBackedBySam(), 0 },
                { basicReadBackedByGoogle(), 0 },
                { new SAMRecordToGATKReadAdapter(samWithISize), 120 },
                { new GoogleGenomicsReadToGATKReadAdapter(googleReadWithFragmentLength), 120 }
        };
    }

    @Test(dataProvider = "GetAndSetFragmentLengthData")
    public void testGetAndSetFragmentLength( final GATKRead read, final int expectedFragmentLength ) {
        Assert.assertEquals(read.getFragmentLength(), expectedFragmentLength, "Wrong fragment length for read");

        read.setFragmentLength(50);
        Assert.assertEquals(read.getFragmentLength(), 50, "Wrong fragment length for read after setFragmentLength()");

        // Negative fragment lengths explicitly allowed
        read.setFragmentLength(-50);
        Assert.assertEquals(read.getFragmentLength(), -50, "Wrong fragment length for read after setFragmentLength()");
    }

    @DataProvider(name = "GetAndSetMappingQualityData")
    public Object[][] getAndSetMappingQualityData() {
        final SAMRecord samWithMappingQualityZero = basicSAMRecord();
        samWithMappingQualityZero.setMappingQuality(0);

        final Read googleReadWithMappingQualityZero = basicGoogleGenomicsRead();
        googleReadWithMappingQualityZero.getAlignment().setMappingQuality(0);

        final Read googleReadWithNoMappingQuality = basicGoogleGenomicsRead();
        googleReadWithNoMappingQuality.getAlignment().setMappingQuality(null);

        return new Object[][]{
                { basicReadBackedBySam(), BASIC_READ_MAPPING_QUALITY },
                { basicReadBackedByGoogle(), BASIC_READ_MAPPING_QUALITY },
                { new SAMRecordToGATKReadAdapter(samWithMappingQualityZero), 0 },
                { new GoogleGenomicsReadToGATKReadAdapter(googleReadWithMappingQualityZero), 0 },
                { new GoogleGenomicsReadToGATKReadAdapter(googleReadWithNoMappingQuality), ReadConstants.NO_MAPPING_QUALITY }
        };
    }

    @Test(dataProvider = "GetAndSetMappingQualityData")
    public void testGetAndSetMappingQuality( final GATKRead read, final int expectedMappingQuality ) {
        Assert.assertEquals(read.getMappingQuality(), expectedMappingQuality, "Wrong mapping quality for read");

        read.setMappingQuality(50);
        Assert.assertEquals(read.getMappingQuality(), 50, "Wrong mapping quality for read after setMappingQuality()");
    }

    @DataProvider(name = "InvalidMappingQualityData")
    public Object[][] invalidMappingQualityData() {
        return new Object[][]{
                { -1 },
                { -10 },
                { 256 },
                { 1000 }
        };
    }

    @Test(dataProvider = "InvalidMappingQualityData", expectedExceptions = IllegalArgumentException.class)
    public void testSetInvalidMappingQualityOnSamBackedRead( final int invalidMappingQuality ) {
        final GATKRead read = basicReadBackedBySam();
        read.setMappingQuality(invalidMappingQuality);
    }

    @Test(dataProvider = "InvalidMappingQualityData", expectedExceptions = IllegalArgumentException.class)
    public void testSetInvalidMappingQualityOnGoogleBackedRead( final int invalidMappingQuality ) {
        final GATKRead read = basicReadBackedByGoogle();
        read.setMappingQuality(invalidMappingQuality);
    }

    @DataProvider(name = "GetAndSetBasesData")
    public Object[][] getAndSetBasesData() {
        final SAMRecord baselessSam = basicSAMRecord();
        baselessSam.setReadBases(null);

        final SAMRecord noAlignedSequenceSam = basicSAMRecord();
        noAlignedSequenceSam.setReadBases(SAMRecord.NULL_SEQUENCE);

        final Read baselessGoogleRead = basicGoogleGenomicsRead();
        baselessGoogleRead.setAlignedSequence(null);

        final Read emptyStringSequenceGoogleRead = basicGoogleGenomicsRead();
        emptyStringSequenceGoogleRead.setAlignedSequence("");

        final Read noAlignedSequenceGoogleRead = basicGoogleGenomicsRead();
        noAlignedSequenceGoogleRead.setAlignedSequence(SAMRecord.NULL_SEQUENCE_STRING);

        return new Object[][]{
                { basicReadBackedBySam(), BASIC_READ_BASES, "ACGT" },
                { basicReadBackedByGoogle(), BASIC_READ_BASES, "ACGT" },
                { new SAMRecordToGATKReadAdapter(baselessSam), new byte[0], "*" },
                { new SAMRecordToGATKReadAdapter(noAlignedSequenceSam), new byte[0], "*" },
                { new GoogleGenomicsReadToGATKReadAdapter(baselessGoogleRead), new byte[0], "*" },
                { new GoogleGenomicsReadToGATKReadAdapter(emptyStringSequenceGoogleRead), new byte[0], "*" },
                { new GoogleGenomicsReadToGATKReadAdapter(noAlignedSequenceGoogleRead), new byte[0], "*" }
        };
    }

    @Test(dataProvider = "GetAndSetBasesData")
    public void testGetAndSetBases( final GATKRead read, final byte[] expectedBases, final String expectedBasesString ) {
        Assert.assertEquals(read.getBases(), expectedBases, "Wrong bases for read");
        Assert.assertEquals(read.getBasesString(), expectedBasesString, "Wrong base string for read");

        final byte[] newBases = {'G', 'C', 'G', 'G'};
        read.setBases(newBases);
        Assert.assertEquals(read.getBases(), newBases, "Wrong bases for read after setBases()");
        Assert.assertEquals(read.getBasesString(), "GCGG", "Wrong base string for read after setBases()");
    }

    @DataProvider(name = "GetAndSetBaseQualitiesData")
    public Object[][] getAndSetBaseQualitiesData() {
        final SAMRecord noQualsSam = basicSAMRecord();
        noQualsSam.setBaseQualities(null);

        final SAMRecord emptyQualsSam = basicSAMRecord();
        emptyQualsSam.setBaseQualities(new byte[0]);

        final Read noQualsGoogleRead = basicGoogleGenomicsRead();
        noQualsGoogleRead.setAlignedQuality(null);

        final Read emptyQualsGoogleRead = basicGoogleGenomicsRead();
        emptyQualsGoogleRead.setAlignedQuality(new ArrayList<>());

        return new Object[][]{
                { basicReadBackedBySam(), BASIC_READ_BASE_QUALITIES },
                { basicReadBackedByGoogle(), BASIC_READ_BASE_QUALITIES },
                { new SAMRecordToGATKReadAdapter(noQualsSam), new byte[0] },
                { new SAMRecordToGATKReadAdapter(emptyQualsSam), new byte[0] },
                { new GoogleGenomicsReadToGATKReadAdapter(noQualsGoogleRead), new byte[0] },
                { new GoogleGenomicsReadToGATKReadAdapter(emptyQualsGoogleRead), new byte[0] }
        };
    }

    @Test(dataProvider = "GetAndSetBaseQualitiesData")
    public void testGetAndSetBaseQualities( final GATKRead read, final byte[] expectedQuals ) {
        Assert.assertEquals(read.getBaseQualities(), expectedQuals, "Wrong base qualities for read");

        final byte[] newQuals = {1, 2, 3, 4};
        read.setBaseQualities(newQuals);
        Assert.assertEquals(read.getBaseQualities(), newQuals, "Wrong base qualities for read after setBaseQualities()");
    }

    @Test(expectedExceptions = GATKException.class)
    public void testSetInvalidBaseQualitiesOnGoogleRead() {
        final GATKRead read = basicReadBackedByGoogle();
        read.setBaseQualities(new byte[]{-1});
    }

    @DataProvider(name = "GetAndSetCigarData")
    public Object[][] getAndSetCigarData() {

        SAMRecord noCigarSam = basicSAMRecord();
        noCigarSam.setCigar(null);

        SAMRecord emptyCigarSam = basicSAMRecord();
        emptyCigarSam.setCigar(new Cigar());

        Read noCigarRead = basicGoogleGenomicsRead();
        noCigarRead.getAlignment().setCigar(null);

        Read emptyCigarRead = basicGoogleGenomicsRead();
        emptyCigarRead.getAlignment().setCigar(null);

        return new Object[][]{
                { basicReadBackedBySam(), TextCigarCodec.decode(BASIC_READ_CIGAR) },
                { basicReadBackedByGoogle(), TextCigarCodec.decode(BASIC_READ_CIGAR) },
                { new SAMRecordToGATKReadAdapter(noCigarSam), new Cigar() },
                { new SAMRecordToGATKReadAdapter(emptyCigarSam), new Cigar() },
                { new GoogleGenomicsReadToGATKReadAdapter(noCigarRead), new Cigar() },
                { new GoogleGenomicsReadToGATKReadAdapter(emptyCigarRead), new Cigar() }
        };
    }

    @Test(dataProvider = "GetAndSetCigarData")
    public void testGetAndSetCigar( final GATKRead read, final Cigar expectedCigar ) {
        Assert.assertEquals(read.getCigar(), expectedCigar, "Wrong cigar for read");

        final Cigar newCigar = TextCigarCodec.decode("4M");
        read.setCigar(newCigar);
        Assert.assertEquals(read.getCigar(), newCigar, "Wrong cigar for read after setCigar()");

        read.setCigar("2M2I");
        Assert.assertEquals(read.getCigar(), TextCigarCodec.decode("2M2I"), "Wrong cigar for read after setCigar()");

        read.setCigar(new Cigar());
        Assert.assertEquals(read.getCigar(), new Cigar(), "Wrong cigar for read after setCigar()");

        read.setCigar((Cigar)null);
        Assert.assertEquals(read.getCigar(), new Cigar(), "Wrong cigar for read after setCigar()");
    }

    @DataProvider(name = "GetAndSetReadGroupData")
    public Object[][] getAndSetReadGroupData() {
        SAMRecord noRGSam = basicSAMRecord();
        noRGSam.clearAttributes();

        Read noRGGoogleRead = basicGoogleGenomicsRead();
        noRGGoogleRead.setReadGroupId(null);

        return new Object[][] {
                { basicReadBackedBySam(), BASIC_READ_GROUP },
                { basicReadBackedByGoogle(), BASIC_READ_GROUP },
                { new SAMRecordToGATKReadAdapter(noRGSam), null },
                { new GoogleGenomicsReadToGATKReadAdapter(noRGGoogleRead), null }
        };
    }

    @Test(dataProvider = "GetAndSetReadGroupData")
    public void testGetAndSetReadGroup( final GATKRead read, final String expectedReadGroup ) {
        Assert.assertEquals(read.getReadGroup(), expectedReadGroup, "Wrong read group for read");

        read.setReadGroup("NewReadGroup");
        Assert.assertEquals(read.getReadGroup(), "NewReadGroup", "Wrong read group for read after setReadGroup()");

        read.setReadGroup(null);
        Assert.assertEquals(read.getReadGroup(), null, "Read group should be null after setReadGroup(null)");
    }

    @DataProvider(name = "IsPairedData")
    public Object[][] isPairedData() {
        SAMRecord unpairedSAM = basicSAMRecord();
        unpairedSAM.setReadPairedFlag(false);

        Read unpairedGoogleRead = basicGoogleGenomicsRead();
        unpairedGoogleRead.setNumberReads(1);

        SAMRecord properlyPairedSAM = basicSAMRecord();
        properlyPairedSAM.setProperPairFlag(true);

        Read properlyPairedGoogleRead = basicGoogleGenomicsRead();
        properlyPairedGoogleRead.setProperPlacement(true);

        SAMRecord unpairedProperlyPairedSAM = basicSAMRecord();
        unpairedProperlyPairedSAM.setReadPairedFlag(false);
        unpairedProperlyPairedSAM.setProperPairFlag(true);

        Read unpairedProperlyPairedGoogleRead = basicGoogleGenomicsRead();
        unpairedProperlyPairedGoogleRead.setNumberReads(1);
        unpairedProperlyPairedGoogleRead.setProperPlacement(true);

        return new Object[][] {
                { basicReadBackedBySam(), true, false },
                { basicReadBackedByGoogle(), true, false },
                { new SAMRecordToGATKReadAdapter(unpairedSAM), false, false },
                { new GoogleGenomicsReadToGATKReadAdapter(unpairedGoogleRead), false, false },
                { new SAMRecordToGATKReadAdapter(properlyPairedSAM), true, true },
                { new GoogleGenomicsReadToGATKReadAdapter(properlyPairedGoogleRead), true, true },

                // We only consider reads to be properly paired if they are also marked as paired
                { new SAMRecordToGATKReadAdapter(unpairedProperlyPairedSAM), false, false },
                { new GoogleGenomicsReadToGATKReadAdapter(unpairedProperlyPairedGoogleRead), false, false }
        };
    }

    @Test(dataProvider = "IsPairedData")
    public void testIsPaired( final GATKRead read, final boolean expectedIsPaired, final boolean expectedIsProperlyPaired ) {
        Assert.assertEquals(read.isPaired(), expectedIsPaired, "Read paired status incorrect");
        Assert.assertEquals(read.isProperlyPaired(), expectedIsProperlyPaired, "Read properly paired status incorrect");

        read.setIsPaired(false);
        Assert.assertFalse(read.isPaired(), "Read should be unpaired after setIsPaired(false)");
        Assert.assertFalse(read.isProperlyPaired(), "Read should not be marked as properly paired after setIsPaired(false)");

        read.setIsPaired(true);
        Assert.assertTrue(read.isPaired(), "Read should be paired after setIsPaired(true)");
        Assert.assertFalse(read.isProperlyPaired(), "Read incorrectly marked as properly paired");

        read.setIsProperlyPaired(true);
        Assert.assertTrue(read.isProperlyPaired(), "Read should be marked as properly paired after setIsProperlyPaired(true)");
        Assert.assertTrue(read.isPaired(), "Read should still be paired after setIsProperlyPaired(true)");

        read.setIsProperlyPaired(false);
        Assert.assertFalse(read.isProperlyPaired(), "Read should not be marked as properly paired after setIsProperlyPaired(false)");
        Assert.assertTrue(read.isPaired(), "Read should still be paired after setIsProperlyPaired(false)");

        read.setIsPaired(false);
        read.setIsProperlyPaired(true);
        Assert.assertTrue(read.isPaired(), "setIsProperlyPaired(true) should have also marked the read as paired");
    }

    @DataProvider(name = "IsUnmappedData")
    public Object[][] isUnmappedData() {
        List<Object[]> testCases = new ArrayList<>();

        for ( final GATKRead unmappedRead : getUnmappedReads() ) {
            testCases.add(new Object[]{ unmappedRead, true });
        }

        testCases.add(new Object[]{ basicReadBackedBySam(), false });
        testCases.add(new Object[]{ basicReadBackedByGoogle(), false });

        return testCases.toArray(new Object[][]{});
    }

    @Test(dataProvider = "IsUnmappedData")
    public void testIsUnmapped( final GATKRead read, final boolean expectedIsUnmapped ) {
        Assert.assertEquals(read.isUnmapped(), expectedIsUnmapped, "Read mapped/unmapped status reported by isUnmapped() is incorrect");

        read.setIsUnmapped();
        Assert.assertTrue(read.isUnmapped(), "Read should be unmapped after setIsUnmapped()");

        read.setPosition("1", 1);
        Assert.assertFalse(read.isUnmapped(), "Read should be mapped after setPosition()");

        read.setIsUnmapped();
        Assert.assertTrue(read.isUnmapped(), "Read should be unmapped after setIsUnmapped()");
    }

    @DataProvider(name = "MateIsUnmappedData")
    public Object[][] mateIsUnmappedData() {
        SAMRecord samWithUnmappedMate = basicSAMRecord();
        samWithUnmappedMate.setMateUnmappedFlag(true);

        SAMRecord samWithUnmappedMate2 = basicSAMRecord();
        samWithUnmappedMate2.setMateReferenceName(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME);

        SAMRecord samWithUnmappedMate3 = basicSAMRecord();
        samWithUnmappedMate3.setMateAlignmentStart(SAMRecord.NO_ALIGNMENT_START);

        Read googleReadWithUnmappedMate = basicGoogleGenomicsRead();
        googleReadWithUnmappedMate.setNextMatePosition(null);

        Read googleReadWithUnmappedMate2 = basicGoogleGenomicsRead();
        googleReadWithUnmappedMate2.getNextMatePosition().setReferenceName(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME);

        Read googleReadWithUnmappedMate3 = basicGoogleGenomicsRead();
        googleReadWithUnmappedMate3.getNextMatePosition().setPosition(-1l);

        Read googleReadWithUnmappedMate4 = basicGoogleGenomicsRead();
        googleReadWithUnmappedMate4.getNextMatePosition().setReferenceName(null);

        Read googleReadWithUnmappedMate5 = basicGoogleGenomicsRead();
        googleReadWithUnmappedMate5.getNextMatePosition().setPosition(null);

        return new Object[][] {
                { basicReadBackedBySam(), false },
                { basicReadBackedByGoogle(), false },
                { new SAMRecordToGATKReadAdapter(samWithUnmappedMate), true },
                { new SAMRecordToGATKReadAdapter(samWithUnmappedMate2), true },
                { new SAMRecordToGATKReadAdapter(samWithUnmappedMate3), true },
                { new GoogleGenomicsReadToGATKReadAdapter(googleReadWithUnmappedMate), true },
                { new GoogleGenomicsReadToGATKReadAdapter(googleReadWithUnmappedMate2), true },
                { new GoogleGenomicsReadToGATKReadAdapter(googleReadWithUnmappedMate3), true },
                { new GoogleGenomicsReadToGATKReadAdapter(googleReadWithUnmappedMate4), true },
                { new GoogleGenomicsReadToGATKReadAdapter(googleReadWithUnmappedMate5), true }
        };
    }

    @Test(dataProvider = "MateIsUnmappedData")
    public void testMateIsUnmapped( final GATKRead read, final boolean expectedMateIsUnmapped ) {
        Assert.assertEquals(read.mateIsUnmapped(), expectedMateIsUnmapped, "Mate unmapped status incorrect");

        read.setMateIsUnmapped();
        Assert.assertTrue(read.mateIsUnmapped(), "Mate should be unmapped after setMateIsUnmapped()");

        read.setMatePosition("1", 1);
        Assert.assertFalse(read.mateIsUnmapped(), "Mate should be mapped after setMatePosition()");

        read.setMateIsUnmapped();
        Assert.assertTrue(read.mateIsUnmapped(), "Mate should be unmapped after setMateIsUnmapped()");

        // Calling setMateIsUnmapped() should have the side effect of marking the read as paired
        read.setIsPaired(false);
        read.setMateIsUnmapped();
        Assert.assertTrue(read.isPaired(), "Read should be marked as paired after call to setMateIsUnmapped()");
    }

    @DataProvider(name = "InvalidMateIsUnmappedData")
    public Object[][] invalidMateIsUnmappedData() {
        SAMRecord unpairedSAM = basicSAMRecord();
        unpairedSAM.setReadPairedFlag(false);

        Read unpairedGoogleRead = basicGoogleGenomicsRead();
        unpairedGoogleRead.setNumberReads(1);

        return new Object[][] {
                { new SAMRecordToGATKReadAdapter(unpairedSAM) },
                { new GoogleGenomicsReadToGATKReadAdapter(unpairedGoogleRead) }
        };
    }

    @Test(dataProvider = "InvalidMateIsUnmappedData", expectedExceptions = IllegalStateException.class)
    public void testInvalidMateIsUnmapped( final GATKRead read ) {
        // Calling mateIsUnmapped() on an unpaired read should throw an error
        read.mateIsUnmapped();
    }

    @DataProvider(name = "IsReverseStrandData")
    public Object[][] isReverseStrandData() {
        SAMRecord reverseStrandSam = basicSAMRecord();
        reverseStrandSam.setReadNegativeStrandFlag(true);

        Read reverseStrandGoogleRead = basicGoogleGenomicsRead();
        reverseStrandGoogleRead.getAlignment().getPosition().setReverseStrand(true);

        return new Object[][] {
                { basicReadBackedBySam(), false },
                { basicReadBackedByGoogle(), false },
                { new SAMRecordToGATKReadAdapter(reverseStrandSam), true },
                { new GoogleGenomicsReadToGATKReadAdapter(reverseStrandGoogleRead), true }
        };
    }

    @Test(dataProvider = "IsReverseStrandData")
    public void testIsReverseStrand( final GATKRead read, final boolean expectedIsReverseStrand ) {
        Assert.assertEquals(read.isReverseStrand(), expectedIsReverseStrand, "Read strand incorrect");

        read.setIsReverseStrand(true);
        Assert.assertTrue(read.isReverseStrand(), "Read should be marked as reverse strand after setIsReverseStrand(true)");

        read.setIsReverseStrand(false);
        Assert.assertFalse(read.isReverseStrand(), "Read should not be marked as reverse strand after setIsReverseStrand(false)");
    }

    @DataProvider(name = "MateIsReverseStrandData")
    public Object[][] mateIsReverseStrandData() {
        SAMRecord samWithReverseStrandMate = basicSAMRecord();
        samWithReverseStrandMate.setMateNegativeStrandFlag(true);

        Read googleReadWithReverseStrandMate = basicGoogleGenomicsRead();
        googleReadWithReverseStrandMate.getNextMatePosition().setReverseStrand(true);

        return new Object[][] {
                { basicReadBackedBySam(), false },
                { basicReadBackedByGoogle(), false },
                { new SAMRecordToGATKReadAdapter(samWithReverseStrandMate), true },
                { new GoogleGenomicsReadToGATKReadAdapter(googleReadWithReverseStrandMate), true }
        };
    }

    @Test(dataProvider = "MateIsReverseStrandData")
    public void testMateIsReverseStrand( final GATKRead read, final boolean expectedMateIsReverseStrand ) {
        Assert.assertEquals(read.mateIsReverseStrand(), expectedMateIsReverseStrand, "Mate strandedness incorrect");

        read.setMateIsReverseStrand(true);
        Assert.assertTrue(read.mateIsReverseStrand(), "Mate should be reverse strand after setMateIsReverseStrand(true)");

        read.setMateIsReverseStrand(false);
        Assert.assertFalse(read.mateIsReverseStrand(), "Mate should not be reverse strand after setMateIsReverseStrand(false)");

        // Calling setMateIsReverseStrand() should have the side effect of marking the read as paired.
        read.setIsPaired(false);
        read.setMateIsReverseStrand(true);
        Assert.assertTrue(read.isPaired(), "Read should be marked as paired after call to setMateIsReverseStrand()");
    }

    @DataProvider(name = "InvalidMateIsReverseStrandData")
    public Object[][] invalidMateIsReverseStrandData() {
        SAMRecord unpairedSAM = basicSAMRecord();
        unpairedSAM.setReadPairedFlag(false);

        Read unpairedGoogleRead = basicGoogleGenomicsRead();
        unpairedGoogleRead.setNumberReads(1);

        return new Object[][] {
                { new SAMRecordToGATKReadAdapter(unpairedSAM) },
                { new GoogleGenomicsReadToGATKReadAdapter(unpairedGoogleRead) }
        };
    }

    @Test(dataProvider = "InvalidMateIsReverseStrandData", expectedExceptions = IllegalStateException.class)
    public void testInvalidMateIsReverseStrand( final GATKRead read ) {
        // Should throw when called on an unpaired read
        read.mateIsReverseStrand();
    }

    @DataProvider(name = "ReadNumberTestData")
    public Object[][] readNumberTestData() {
        SAMRecord secondOfPairSam = basicSAMRecord();
        secondOfPairSam.setSecondOfPairFlag(true);
        secondOfPairSam.setFirstOfPairFlag(false);

        Read secondOfPairGoogleRead = basicGoogleGenomicsRead();
        secondOfPairGoogleRead.setReadNumber(1);

        return new Object[][] {
                { basicReadBackedBySam(), true, false },
                { basicReadBackedByGoogle(), true, false },
                { new SAMRecordToGATKReadAdapter(secondOfPairSam), false, true },
                { new GoogleGenomicsReadToGATKReadAdapter(secondOfPairGoogleRead), false, true }
        };
    }

    @Test(dataProvider = "ReadNumberTestData")
    public void testReadNumber( final GATKRead read, final boolean expectedIsFirstOfPair, final boolean expectedIsSecondOfPair ) {
        Assert.assertEquals(read.isFirstOfPair(), expectedIsFirstOfPair, "Read first of pair setting incorrect");
        Assert.assertEquals(read.isSecondOfPair(), expectedIsSecondOfPair, "Read second of pair setting incorrect");

        read.setIsFirstOfPair();
        Assert.assertTrue(read.isFirstOfPair(), "Read should be marked first of pair after setIsFirstOfPair()");
        Assert.assertFalse(read.isSecondOfPair(), "Read should not be marked second of pair after setIsFirstOfPair()");

        read.setIsSecondOfPair();
        Assert.assertFalse(read.isFirstOfPair(), "Read should not be marked first of pair after setIsSecondOfPair()");
        Assert.assertTrue(read.isSecondOfPair(), "Read should be marked second of pair after setIsSecondOfPair()");
    }

    @DataProvider(name = "GetReadNumberOfUnpairedReadData")
    public Object[][] getReadNumberOfUnpairedReadData() {
        SAMRecord unpairedSAM = basicSAMRecord();
        unpairedSAM.setFirstOfPairFlag(true);
        unpairedSAM.setSecondOfPairFlag(false);
        unpairedSAM.setReadPairedFlag(false);

        SAMRecord unpairedSAM2 = basicSAMRecord();
        unpairedSAM2.setSecondOfPairFlag(true);
        unpairedSAM2.setFirstOfPairFlag(false);
        unpairedSAM2.setReadPairedFlag(false);

        Read unpairedGoogleRead = basicGoogleGenomicsRead();
        unpairedGoogleRead.setReadNumber(0);
        unpairedGoogleRead.setNumberReads(1);

        Read unpairedGoogleRead2 = basicGoogleGenomicsRead();
        unpairedGoogleRead2.setReadNumber(1);
        unpairedGoogleRead2.setNumberReads(1);

        return new Object[][] {
                { new SAMRecordToGATKReadAdapter(unpairedSAM) },
                { new SAMRecordToGATKReadAdapter(unpairedSAM2) },
                { new GoogleGenomicsReadToGATKReadAdapter(unpairedGoogleRead) },
                { new GoogleGenomicsReadToGATKReadAdapter(unpairedGoogleRead2) }
        };
    }

    @Test(dataProvider = "GetReadNumberOfUnpairedReadData")
    public void testGetReadNumberOfUnpairedRead( final GATKRead unpairedRead ) {
        Assert.assertFalse(unpairedRead.isFirstOfPair(), "Unpaired read should not be marked as first of pair");
        Assert.assertFalse(unpairedRead.isSecondOfPair(), "Unpaired read should not be marked as second of pair");
    }

    @DataProvider(name = "GetAndSetSimpleFlagsData")
    public Object[][] getAndSetSimpleFlagsData() {
        return new Object[][] {
                { basicReadBackedBySam() },
                { basicReadBackedByGoogle() }
        };
    }

    @Test(dataProvider = "GetAndSetSimpleFlagsData")
    public void testGetAndSetSimpleFlags( final GATKRead read ) {
        read.setIsSecondaryAlignment(true);
        Assert.assertTrue(read.isSecondaryAlignment(), "Read should be marked as secondary alignment after setIsSecondaryAlignment(true)");

        read.setIsSecondaryAlignment(false);
        Assert.assertFalse(read.isSecondaryAlignment(), "Read should not be marked as secondary alignment after setIsSecondaryAlignment(false)");

        read.setIsSupplementaryAlignment(true);
        Assert.assertTrue(read.isSupplementaryAlignment(), "Read should be marked as supplementary alignment after setIsSupplementaryAlignment(true)");

        read.setIsSupplementaryAlignment(false);
        Assert.assertFalse(read.isSupplementaryAlignment(), "Read should not be marked as supplementary alignment after setIsSupplementaryAlignment(false)");

        read.setFailsVendorQualityCheck(true);
        Assert.assertTrue(read.failsVendorQualityCheck(), "Read should be marked as failing vendor quality checks after setFailsVendorQualityCheck(true)");

        read.setFailsVendorQualityCheck(false);
        Assert.assertFalse(read.failsVendorQualityCheck(), "Read should not be marked as failing vendor quality checks after setFailsVendorQualityCheck(false)");

        read.setIsDuplicate(true);
        Assert.assertTrue(read.isDuplicate(), "Read should be marked as a duplicate after setIsDuplicate(true)");

        read.setIsDuplicate(false);
        Assert.assertFalse(read.isDuplicate(), "Read should not be marked as a duplicate after setIsDuplicate(false)");
    }

    @DataProvider(name = "GetAndSetAttributesData")
    public Object[][] getAndSetAttributesData() {
        return new Object[][] {
                { basicReadBackedBySam() },
                { basicReadBackedByGoogle() }
        };
    }

    @Test(dataProvider = "GetAndSetAttributesData")
    public void testGetAndSetIntegerAttribute( final GATKRead read ) {
        Assert.assertNull(read.getAttributeAsInteger("DR"), "Attribute DR should be null");

        read.setAttribute("DR", 5);
        Assert.assertEquals(read.getAttributeAsInteger("DR").intValue(), 5, "Wrong value for attribute DR");

        // test type coercion
        read.setAttribute("DR", "6");
        Assert.assertEquals(read.getAttributeAsInteger("DR").intValue(), 6, "Wrong value for attribute DR");

        read.setAttribute("DR", 10);
        Assert.assertEquals(read.getAttributeAsInteger("DR").intValue(), 10, "Wrong value for attribute DR");

        read.clearAttribute("DR");
        Assert.assertNull(read.getAttributeAsInteger("DR"), "Attribute DR should be null");
    }

    @Test(dataProvider = "GetAndSetAttributesData", expectedExceptions = GATKException.ReadAttributeTypeMismatch.class)
    public void testGetNonIntegerAttributeAsInteger( final GATKRead read ) {
        read.setAttribute("DR", "notaninteger");
        read.getAttributeAsInteger("DR");
    }

    @Test(dataProvider = "GetAndSetAttributesData")
    public void testGetAndSetStringAttribute( final GATKRead read ) {
        Assert.assertNull(read.getAttributeAsString("DR"), "Attribute DR should be null");

        read.setAttribute("DR", "foo");
        Assert.assertEquals(read.getAttributeAsString("DR"), "foo", "Wrong value for attribute DR");

        // Test type coercion
        read.setAttribute("DR", 5);
        Assert.assertEquals(read.getAttributeAsString("DR"), "5", "Wrong value for attribute DR");

        read.setAttribute("DR", "bar");
        Assert.assertEquals(read.getAttributeAsString("DR"), "bar", "Wrong value for attribute DR");

        read.clearAttribute("DR");
        Assert.assertNull(read.getAttributeAsString("DR"), "Attribute DR should be null");
    }

    @Test(dataProvider = "GetAndSetAttributesData")
    public void testGetAndSetByteArrayAttribute( final GATKRead read ) {
        Assert.assertNull(read.getAttributeAsByteArray("DR"), "Attribute DR should be null");

        read.setAttribute("DR", new byte[]{1, 2, 3});
        Assert.assertEquals(read.getAttributeAsByteArray("DR"), new byte[]{1, 2, 3}, "Wrong value for attribute DR");

        // Test type coercion
        read.setAttribute("DR", "abc");
        Assert.assertEquals(read.getAttributeAsByteArray("DR"), new byte[]{'a', 'b', 'c'}, "Wrong value for attribute DR");

        read.setAttribute("DR", new byte[]{4, 5, 6, 7});
        Assert.assertEquals(read.getAttributeAsByteArray("DR"), new byte[]{4, 5, 6, 7}, "Wrong value for attribute DR");

        read.clearAttribute("DR");
        Assert.assertNull(read.getAttributeAsByteArray("DR"), "Attribute DR should be null");
    }

    @Test(expectedExceptions = GATKException.ReadAttributeTypeMismatch.class)
    public void testGetNonByteArrayAttributeAsByteArray() {
        // Only SAMRecord-backed reads can encounter a type mismatch for byte array attributes
        SAMRecord sam = basicSAMRecord();
        sam.setAttribute("DR", 5);
        final GATKRead samBackedRead = new SAMRecordToGATKReadAdapter(sam);

        samBackedRead.getAttributeAsByteArray("DR");
    }

    @Test(dataProvider = "GetAndSetAttributesData")
    public void testClearAllAttributes( final GATKRead read ) {
        read.setAttribute("DR", "foo");
        read.setAttribute("LT", 5);

        Assert.assertNotNull(read.getAttributeAsString("DR"), "Attribute DR should not be null");
        Assert.assertNotNull(read.getAttributeAsInteger("LT"), "Attribute LT should not be null");

        read.clearAttributes();
        Assert.assertNull(read.getAttributeAsString("DR"), "Attribute DR should be null");
        Assert.assertNull(read.getAttributeAsInteger("LT"), "Attribute LT should be null");
    }
}