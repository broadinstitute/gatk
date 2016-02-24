package org.broadinstitute.hellbender.utils.pileup;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

import static org.broadinstitute.hellbender.utils.read.ReadConstants.NO_MAPPING_QUALITY;

/**
 * Test routines for read-backed pileup.
 */
public final class ReadPileupUnitTest {
    private static SAMFileHeader header;
    private Locatable loc;

    @BeforeClass
    public void beforeClass() {
        header = ArtificialReadUtils.createArtificialSamHeader(1, 1, 1000);
        loc = new SimpleInterval("chr1", 1, 1);
    }

    /**
     * Ensure that basic read group splitting works.
     */
    @Test
    public void testSplitByReadGroup() {
        final SAMReadGroupRecord readGroupOne = new SAMReadGroupRecord("rg1");
        final SAMReadGroupRecord readGroupTwo = new SAMReadGroupRecord("rg2");

        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader(1, 1, 1000);
        header.addReadGroup(readGroupOne);
        header.addReadGroup(readGroupTwo);

        final GATKRead read1 = ArtificialReadUtils.createArtificialRead(header,"read1",0,1,10);
        read1.setReadGroup(readGroupOne.getId());
        final GATKRead read2 = ArtificialReadUtils.createArtificialRead(header,"read2",0,1,10);
        read2.setReadGroup(readGroupTwo.getId());
        final GATKRead read3 = ArtificialReadUtils.createArtificialRead(header,"read3",0,1,10);
        read3.setReadGroup(readGroupOne.getId());
        final GATKRead read4 = ArtificialReadUtils.createArtificialRead(header,"read4",0,1,10);
        read4.setReadGroup(readGroupTwo.getId());
        final GATKRead read5 = ArtificialReadUtils.createArtificialRead(header,"read5",0,1,10);
        read5.setReadGroup(readGroupTwo.getId());
        final GATKRead read6 = ArtificialReadUtils.createArtificialRead(header,"read6",0,1,10);
        read6.setReadGroup(readGroupOne.getId());
        final GATKRead read7 = ArtificialReadUtils.createArtificialRead(header,"read7",0,1,10);
        read7.setReadGroup(readGroupOne.getId());

        final ReadPileup pileup = new ReadPileup(loc, Arrays.asList(read1, read2, read3, read4, read5, read6, read7), 1);

        final ReadPileup rg1Pileup = pileup.makeFilteredPileup(pe -> "rg1".equals(pe.getRead().getReadGroup()));
        final List<GATKRead> rg1Reads = rg1Pileup.getReads();
        Assert.assertEquals(rg1Reads.size(), 4, "Wrong number of reads in read group rg1");
        Assert.assertEquals(rg1Reads.get(0), read1, "Read " + read1.getName() + " should be in rg1 but isn't");
        Assert.assertEquals(rg1Reads.get(1), read3, "Read " + read3.getName() + " should be in rg1 but isn't");
        Assert.assertEquals(rg1Reads.get(2), read6, "Read " + read6.getName() + " should be in rg1 but isn't");
        Assert.assertEquals(rg1Reads.get(3), read7, "Read " + read7.getName() + " should be in rg1 but isn't");

        final ReadPileup rg2Pileup = pileup.makeFilteredPileup(pe -> "rg2".equals(pe.getRead().getReadGroup()));
        final List<GATKRead> rg2Reads = rg2Pileup.getReads();
        Assert.assertEquals(rg2Reads.size(), 3, "Wrong number of reads in read group rg2");
        Assert.assertEquals(rg2Reads.get(0), read2, "Read " + read2.getName() + " should be in rg2 but isn't");
        Assert.assertEquals(rg2Reads.get(1), read4, "Read " + read4.getName() + " should be in rg2 but isn't");
        Assert.assertEquals(rg2Reads.get(2), read5, "Read " + read5.getName() + " should be in rg2 but isn't");
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testNullConstructorParametersReads() {
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader(1, 1, 1000);

        ArtificialReadUtils.createArtificialRead(header,"read1",0,1,10);

        new ReadPileup(loc, null, Arrays.asList(1));
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testNullConstructorParametersOffsets() {
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader(1, 1, 1000);

        final GATKRead read1 = ArtificialReadUtils.createArtificialRead(header, "read1", 0, 1, 10);

        new ReadPileup(loc, Arrays.asList(read1), null);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testNullConstructorParametersLoc() {
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader(1, 1, 1000);

        final GATKRead read1 = ArtificialReadUtils.createArtificialRead(header,"read1",0,1,10);

        new ReadPileup(null, Arrays.asList(read1), Arrays.asList(1));
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testInvalidConstructorParametersReadsAndOffsetsLists() {
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader(1, 1, 1000);

        final GATKRead read1 = ArtificialReadUtils.createArtificialRead(header, "read1", 0, 1, 10);

        new ReadPileup(loc, Arrays.asList(read1), Arrays.asList(1, 2, 3));
    }

    /**
     * Ensure that splitting read groups still works when dealing with null read groups.
     */
    @Test
    public void testSplitByNullReadGroups() {
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader(1,1,1000);

        final GATKRead read1 = ArtificialReadUtils.createArtificialRead(header,"read1",0,1,10);
        final GATKRead read2 = ArtificialReadUtils.createArtificialRead(header,"read2",0,1,10);
        final GATKRead read3 = ArtificialReadUtils.createArtificialRead(header,"read3",0,1,10);

        final ReadPileup pileup = new ReadPileup(loc,
                                                           Arrays.asList(read1, read2, read3),
                                                           Arrays.asList(1, 1, 1));

        final ReadPileup nullRgPileup = pileup.makeFilteredPileup(pe -> pe.getRead().getReadGroup() == null);
        final List<GATKRead> nullRgReads = nullRgPileup.getReads();
        Assert.assertEquals(nullRgPileup.size(), 3, "Wrong number of reads in null read group");
        Assert.assertEquals(nullRgReads.get(0), read1, "Read " + read1.getName() + " should be in null rg but isn't");
        Assert.assertEquals(nullRgReads.get(1), read2, "Read " + read2.getName() + " should be in null rg but isn't");
        Assert.assertEquals(nullRgReads.get(2), read3, "Read " + read3.getName() + " should be in null rg but isn't");

        final ReadPileup rg1Pileup = pileup.makeFilteredPileup(pe -> "rg1".equals(pe.getRead().getReadGroup()));
        Assert.assertTrue(rg1Pileup.isEmpty(), "Pileup for non-existent read group should return empty pileup");
    }

    private final class ReadPileupCountTest {
        final String sample;
        final int nReads, nMapq0, nDeletions;

        private ReadPileupCountTest(final int nReads, final int nMapq0, final int nDeletions) {
            this.sample = "sample" + UUID.randomUUID();
            this.nReads = nReads;
            this.nMapq0 = nMapq0;
            this.nDeletions = nDeletions;
        }

        private List<PileupElement> makeReads( final int n, final int mapq, final String op ) {
            final int readLength = 3;

            final List<PileupElement> elts = new ArrayList<>(n);
            for ( int i = 0; i < n; i++ ) {
                final GATKRead read = ArtificialReadUtils.createArtificialRead(header, "read", 0, 1, readLength);
                read.setBases(Utils.dupBytes((byte) 'A', readLength));
                read.setBaseQualities(Utils.dupBytes((byte) 30, readLength));
                read.setCigar("1M1" + op + "1M");
                read.setMappingQuality(mapq);
                final int baseOffset = op.equals("M") ? 1 : 0;
                final CigarElement cigarElement = read.getCigar().getCigarElement(1);
                elts.add(new PileupElement(read, baseOffset, cigarElement, 1, 0));
            }

            return elts;
        }

        private ReadPileup makePileup() {
            final List<PileupElement> elts = new ArrayList<>(nReads);

            elts.addAll(makeReads(nMapq0, 0, "M"));
            elts.addAll(makeReads(nDeletions, 30, "D"));
            elts.addAll(makeReads(nReads - nMapq0 - nDeletions, 30, "M"));

            return new ReadPileup(loc, elts);
        }

        @Override
        public String toString() {
            return "ReadPileupCountTest{" +
                    "sample='" + sample + '\'' +
                    ", nReads=" + nReads +
                    ", nMapq0=" + nMapq0 +
                    ", nDeletions=" + nDeletions +
                    '}';
        }
    }

    @DataProvider(name = "ReadPileupCountingTest")
    public Object[][] makeReadPileupCountingTest() {
        final List<Object[]> tests = new ArrayList<>();

        for ( final int nMapq : Arrays.asList(0, 10, 20) ) {
            for ( final int nDeletions : Arrays.asList(0, 10, 20) ) {
                for ( final int nReg : Arrays.asList(0, 10, 20) ) {
                    final int total = nMapq + nDeletions + nReg;
                    if ( total > 0 ) {
                        tests.add(new Object[]{new ReadPileupCountTest(total, nMapq, nDeletions)});
                    }
                }
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "ReadPileupCountingTest")
    public void testReadPileupCountingTestSinglePileup(final ReadPileupCountTest params) {
        testReadPileupCounts(params.makePileup(), params);
    }

    private void testReadPileupCounts(final ReadPileup pileup, final ReadPileupCountTest expected) {
        for ( int cycles = 0; cycles < 3; cycles++ ) {
            // multiple cycles to make sure caching is working
            Assert.assertEquals(pileup.size(), expected.nReads);
            Assert.assertEquals(pileup.getNumberOfElements(p -> p.isDeletion()), expected.nDeletions);
            Assert.assertEquals(pileup.getNumberOfElements(p -> true), expected.nReads);
            Assert.assertEquals(pileup.getNumberOfElements(p -> false), 0);
            Assert.assertEquals(pileup.getNumberOfElements(p -> p.getRead().getMappingQuality() == 0), expected.nMapq0);
        }
    }


    @Test
    public void testReadPileupMappingQuals() {

        // create a read with high MQ
        final GATKRead read = ArtificialReadUtils.createArtificialRead(header, "read", 0, 1, 10);
        read.setBases(Utils.dupBytes((byte) 'A', 10));
        read.setBaseQualities(Utils.dupBytes((byte) 30, 10));
        read.setCigar("10M");
        read.setMappingQuality(200); // set a MQ higher than max signed byte

        // now create the pileup
        final List<PileupElement> elts = new LinkedList<>();
        elts.add(new PileupElement(read, 0, read.getCigar().getCigarElement(0), 0, 0));
        final ReadPileup pileup = new ReadPileup(loc, elts);

        Assert.assertEquals(pileup.getMappingQuals()[0], 200);
    }

    @Test
    public void testEmptyPileup(){
        final ReadPileup empty = new ReadPileup(loc);
        Assert.assertTrue(empty.isEmpty());
        Assert.assertTrue(empty.size() == 0);
        Assert.assertEquals(empty.getBaseCounts(), new int[4]);
        Assert.assertEquals(empty.getBases(), new byte[0]);
        Assert.assertEquals(empty.getBaseQuals(), new byte[0]);
        Assert.assertEquals(empty.getLocation(), loc);
        Assert.assertEquals(empty.getMappingQuals(), new int[0]);
        Assert.assertEquals(empty.getNumberOfElements(p -> true), 0);
        Assert.assertEquals(empty.getNumberOfElements(p -> false), 0);
        Assert.assertEquals(empty.getNumberOfElements(p -> p.isDeletion()), 0);
        Assert.assertEquals(empty.getNumberOfElements(p -> p.isBeforeDeletionStart()), 0);
        Assert.assertEquals(empty.getNumberOfElements(p -> p.isBeforeInsertion()), 0);
        Assert.assertEquals(empty.getNumberOfElements(p -> p.getRead().getMappingQuality() == 0), 0);
        Assert.assertTrue(empty.getOffsets().isEmpty());
        Assert.assertTrue(empty.getReadGroupIDs().isEmpty());
        Assert.assertTrue(empty.getReads().isEmpty());
        Assert.assertTrue(empty.getSamples(header).isEmpty());
        Assert.assertTrue(empty.getPileupForLane("fred").isEmpty());
    }

    @Test
    public void testSimplePileup(){
        final int readlength = 10;
        final byte[] bases1 = Utils.repeatChars('A', readlength);
        final byte[] quals1 = Utils.repeatBytes((byte) 10, readlength);
        final String cigar1 = "10M";

        final byte[] bases2 = Utils.repeatChars('C', readlength);
        final byte[] quals2 = Utils.repeatBytes((byte) 20, readlength);
        final String cigar2 = "5M3I2M";

        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader();
        final GATKRead read1 = ArtificialReadUtils.createArtificialRead(header, bases1, quals1, cigar1);
        read1.setName("read1");
        final GATKRead read2 = ArtificialReadUtils.createArtificialRead(header, bases2, quals2, cigar2);
        read1.setName("read2");
        final List<GATKRead> reads = Arrays.asList(read1, read2);
        final ReadPileup pu = new ReadPileup(loc, reads, 1);

        Assert.assertNotNull(pu.toString()); //checking non-blowup. We're not making any claims about the format of toString

        Assert.assertEquals(pu.getBases(), new byte[]{(byte) 'A', (byte) 'C'});
        Assert.assertFalse(pu.isEmpty());
        Assert.assertEquals(pu.size(), 2, "size");
        Assert.assertEquals(pu.getBaseCounts(), new int[]{1, 1, 0, 0});
        Assert.assertEquals(pu.getBaseQuals(), new byte[]{10, 20});
        Assert.assertEquals(pu.getLocation(), loc);
        Assert.assertEquals(pu.getMappingQuals(), new int[]{NO_MAPPING_QUALITY, NO_MAPPING_QUALITY});
        Assert.assertEquals(pu.makeFilteredPileup(p -> p.getQual() >= 12).makeFilteredPileup(p -> p.getMappingQual() >= 0).size(), 1, "getBaseAndMappingFilteredPileup");
        Assert.assertEquals(pu.makeFilteredPileup(r -> r.getQual() >= 12).size(), 1, "getBaseFilteredPileup");
        Assert.assertEquals(pu.getNumberOfElements(p -> true), 2);
        Assert.assertEquals(pu.getNumberOfElements(p -> false), 0);
        Assert.assertEquals(pu.getNumberOfElements(p -> p.isDeletion()), 0);
        Assert.assertEquals(pu.getNumberOfElements(p -> p.isBeforeDeletionStart()), 0);
        Assert.assertEquals(pu.getNumberOfElements(p -> p.isBeforeInsertion()), 0);
        Assert.assertEquals(pu.getNumberOfElements(p -> p.getRead().getMappingQuality() == 0), 2);
        Assert.assertEquals(pu.getOffsets(), Arrays.asList(1, 1), "getOffsets");
        Assert.assertEquals(pu.getReadGroupIDs(), Arrays.asList(new SAMReadGroupRecord[]{null}), "getReadGroups");
        Assert.assertEquals(pu.getReads(), Arrays.asList(read1, read2), "getReads");
        Assert.assertEquals(pu.getSamples(header), Arrays.asList(new String[]{null}), "getSamples");
        Assert.assertTrue(pu.getPileupForLane("fred").isEmpty());
        Assert.assertTrue(pu.makeFilteredPileup(r -> r.getMappingQual() >= 10).isEmpty());
    }

    @Test
    public void testSimplePileupWithIndelsOffset(){
        final int readlength = 10;
        final byte[] bases1 = Utils.repeatChars('A', readlength);
        final byte[] quals1 = Utils.repeatBytes((byte) 10, readlength);
        final String cigar1 = "10M";

        final byte[] bases2 = Utils.repeatChars('C', readlength);
        final byte[] quals2 = Utils.repeatBytes((byte) 20, readlength);
        final String cigar2 = "5M3I2M";

        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader();
        final GATKRead read1 = ArtificialReadUtils.createArtificialRead(header, bases1, quals1, cigar1);
        read1.setName("read1");
        final GATKRead read2 = ArtificialReadUtils.createArtificialRead(header, bases2, quals2, cigar2);
        read1.setName("read2");
        final List<GATKRead> reads = Arrays.asList(read1, read2);

        final int offset = 4;
        final ReadPileup pu = new ReadPileup(loc, reads, offset);

        Assert.assertNotNull(pu.toString()); //checking non-blowup. We're not making any claims about the format of toString

        Assert.assertEquals(pu.getBases(), new byte[]{(byte) 'A', (byte) 'C'});
        Assert.assertFalse(pu.isEmpty());
        Assert.assertEquals(pu.size(), 2, "size");
        Assert.assertEquals(pu.getBaseCounts(), new int[]{1, 1, 0, 0});
        Assert.assertEquals(pu.getBaseQuals(), new byte[]{10, 20});
        Assert.assertEquals(pu.getLocation(), loc);
        Assert.assertEquals(pu.getMappingQuals(), new int[]{NO_MAPPING_QUALITY, NO_MAPPING_QUALITY});
        Assert.assertEquals(pu.makeFilteredPileup(p -> p.getQual() >= 12).makeFilteredPileup(p -> p.getMappingQual() >= 0).size(), 1, "getBaseAndMappingFilteredPileup");
        Assert.assertEquals(pu.makeFilteredPileup(r -> r.getQual() >= 12).size(), 1, "getBaseFilteredPileup");
        Assert.assertEquals(pu.getNumberOfElements(p -> true), 2);
        Assert.assertEquals(pu.getNumberOfElements(p -> false), 0);
        Assert.assertEquals(pu.getNumberOfElements(p -> p.isDeletion()), 0);
        Assert.assertEquals(pu.getNumberOfElements(p -> p.isBeforeDeletionStart()), 0);
        Assert.assertEquals(pu.getNumberOfElements(p -> p.isBeforeInsertion()), 1);
        Assert.assertEquals(pu.getNumberOfElements(p -> p.getRead().getMappingQuality() == 0), 2);
        Assert.assertEquals(pu.getOffsets(), Arrays.asList(offset, offset), "getOffsets");
        Assert.assertEquals(pu.getReadGroupIDs(), Arrays.asList(new SAMReadGroupRecord[]{null}), "getReadGroups");
        Assert.assertEquals(pu.getReads(), Arrays.asList(read1, read2), "getReads");
        Assert.assertEquals(pu.getSamples(header), Arrays.asList(new String[]{null}), "getSamples");
        Assert.assertTrue(pu.getPileupForLane("fred").isEmpty());
        Assert.assertTrue(pu.makeFilteredPileup(r -> r.getMappingQual() >= 10).isEmpty());
    }

    @Test
    public void testSimplePileupWithOffset(){
        final int readlength = 10;
        final byte[] bases1 = Utils.repeatChars('A', readlength);
        final byte[] quals1 = Utils.repeatBytes((byte)10, readlength);
        final String cigar1 = "10M";

        final byte[] bases2 = Utils.repeatChars('C', readlength);
        final byte[] quals2 = Utils.repeatBytes((byte) 20, readlength);
        final String cigar2 = "10M";

        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader();
        final GATKRead read1 = ArtificialReadUtils.createArtificialRead(header, bases1, quals1, cigar1);
        read1.setName("read1");
        final GATKRead read2 = ArtificialReadUtils.createArtificialRead(header, bases2, quals2, cigar2);
        read1.setName("read2");
        final List<GATKRead> reads = Arrays.asList(read1, read2);
        final int off = 6;
        final ReadPileup pu = new ReadPileup(loc, reads, off);
        Assert.assertEquals(pu.getBases(), new byte[]{(byte) 'A', (byte) 'C'});

        Assert.assertNotNull(pu.toString()); //checking non-blowup. We're not making any claims about the format of toString

        Assert.assertFalse(pu.isEmpty());
        Assert.assertEquals(pu.size(), 2, "size");
        Assert.assertEquals(pu.getBaseCounts(), new int[]{1, 1, 0, 0});
        Assert.assertEquals(pu.getBaseQuals(), new byte[]{10, 20});
        Assert.assertEquals(pu.getLocation(), loc);
        Assert.assertEquals(pu.getMappingQuals(), new int[]{NO_MAPPING_QUALITY, NO_MAPPING_QUALITY});
        Assert.assertTrue(pu.makeFilteredPileup(r -> r.getRead().isReverseStrand()).isEmpty(), "getReverseStrandPileup");
        Assert.assertEquals(pu.makeFilteredPileup(r -> !r.getRead().isReverseStrand()).size(), 2, "getForwardStrandPileup");
        Assert.assertEquals(pu.makeFilteredPileup(p -> p.getQual() >= 12).makeFilteredPileup(p -> p.getMappingQual() >= 0).size(), 1, "getBaseAndMappingFilteredPileup");
        Assert.assertEquals(pu.makeFilteredPileup(p -> p.getQual() >= 12).size(), 1, "getBaseFilteredPileup");
        Assert.assertEquals(pu.getNumberOfElements(p -> true), 2);
        Assert.assertEquals(pu.getNumberOfElements(p -> false), 0);
        Assert.assertEquals(pu.getNumberOfElements(p -> p.isDeletion()), 0);
        Assert.assertEquals(pu.getNumberOfElements(p -> p.isBeforeDeletionStart()), 0);
        Assert.assertEquals(pu.getNumberOfElements(p -> p.isBeforeInsertion()), 0);
        Assert.assertEquals(pu.getNumberOfElements(p -> p.getRead().getMappingQuality() == 0), 2);
        Assert.assertEquals(pu.getOffsets(), Arrays.asList(off, off), "getOffsets");
        Assert.assertEquals(pu.getReadGroupIDs(), Arrays.asList(new SAMReadGroupRecord[]{null}), "getReadGroups");
        Assert.assertEquals(pu.getReads(), Arrays.asList(read1, read2), "getReads");
        Assert.assertEquals(pu.getSamples(header), Arrays.asList(new String[]{null}), "getSamples");
        Assert.assertTrue(pu.getPileupForLane("fred").isEmpty());
        Assert.assertTrue(pu.makeFilteredPileup(p -> p.getMappingQual() >= 10).isEmpty());
    }

}