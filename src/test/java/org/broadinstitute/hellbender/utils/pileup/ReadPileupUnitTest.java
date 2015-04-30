package org.broadinstitute.hellbender.utils.pileup;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.ArtificialSAMUtils;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;
import java.util.stream.Collectors;

import static htsjdk.samtools.SAMRecord.NO_MAPPING_QUALITY;

/**
 * Test routines for read-backed pileup.
 */
public class ReadPileupUnitTest {
    protected static SAMFileHeader header;
    private Locatable loc;

    @BeforeClass
    public void beforeClass() {
        header = ArtificialSAMUtils.createArtificialSamHeader(1, 1, 1000);
        loc = new SimpleInterval("chr1", 1, 1);
    }

    /**
     * Ensure that basic read group splitting works.
     */
    @Test
    public void testSplitByReadGroup() {
        SAMReadGroupRecord readGroupOne = new SAMReadGroupRecord("rg1");
        SAMReadGroupRecord readGroupTwo = new SAMReadGroupRecord("rg2");

        SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(1,1,1000);
        header.addReadGroup(readGroupOne);
        header.addReadGroup(readGroupTwo);

        SAMRecord read1 = ArtificialSAMUtils.createArtificialRead(header,"read1",0,1,10);
        read1.setAttribute("RG",readGroupOne.getId());
        SAMRecord read2 = ArtificialSAMUtils.createArtificialRead(header,"read2",0,1,10);
        read2.setAttribute("RG",readGroupTwo.getId());
        SAMRecord read3 = ArtificialSAMUtils.createArtificialRead(header,"read3",0,1,10);
        read3.setAttribute("RG",readGroupOne.getId());
        SAMRecord read4 = ArtificialSAMUtils.createArtificialRead(header,"read4",0,1,10);
        read4.setAttribute("RG",readGroupTwo.getId());
        SAMRecord read5 = ArtificialSAMUtils.createArtificialRead(header,"read5",0,1,10);
        read5.setAttribute("RG",readGroupTwo.getId());
        SAMRecord read6 = ArtificialSAMUtils.createArtificialRead(header,"read6",0,1,10);
        read6.setAttribute("RG",readGroupOne.getId());
        SAMRecord read7 = ArtificialSAMUtils.createArtificialRead(header,"read7",0,1,10);
        read7.setAttribute("RG",readGroupOne.getId());

        ReadPileup pileup = new ReadPileup(null, Arrays.asList(read1, read2, read3, read4, read5, read6, read7), 1);

        ReadPileup rg1Pileup = pileup.getPileupForReadGroup("rg1");
        List<SAMRecord> rg1Reads = rg1Pileup.getReads();
        Assert.assertEquals(rg1Reads.size(), 4, "Wrong number of reads in read group rg1");
        Assert.assertEquals(rg1Reads.get(0), read1, "Read " + read1.getReadName() + " should be in rg1 but isn't");
        Assert.assertEquals(rg1Reads.get(1), read3, "Read " + read3.getReadName() + " should be in rg1 but isn't");
        Assert.assertEquals(rg1Reads.get(2), read6, "Read " + read6.getReadName() + " should be in rg1 but isn't");
        Assert.assertEquals(rg1Reads.get(3), read7, "Read " + read7.getReadName() + " should be in rg1 but isn't");

        ReadPileup rg2Pileup = pileup.getPileupForReadGroup("rg2");
        List<SAMRecord> rg2Reads = rg2Pileup.getReads();
        Assert.assertEquals(rg2Reads.size(), 3, "Wrong number of reads in read group rg2");
        Assert.assertEquals(rg2Reads.get(0), read2, "Read " + read2.getReadName() + " should be in rg2 but isn't");
        Assert.assertEquals(rg2Reads.get(1), read4, "Read " + read4.getReadName() + " should be in rg2 but isn't");
        Assert.assertEquals(rg2Reads.get(2), read5, "Read " + read5.getReadName() + " should be in rg2 but isn't");
    }

    @Test(expectedExceptions = GATKException.class)
    public void testException1() {
        SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(1,1,1000);

        SAMRecord read1 = ArtificialSAMUtils.createArtificialRead(header,"read1",0,1,10);

        ReadPileup pileup = new ReadPileup(null,
                null,
                Arrays.asList(1));
    }

    @Test(expectedExceptions = GATKException.class)
    public void testException2() {
        SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(1, 1, 1000);

        SAMRecord read1 = ArtificialSAMUtils.createArtificialRead(header,"read1",0,1,10);

        ReadPileup pileup = new ReadPileup(null,
                Arrays.asList(read1),
                null);
    }

    @Test(expectedExceptions = GATKException.class)
    public void testException3() {
        SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(1, 1, 1000);

        SAMRecord read1 = ArtificialSAMUtils.createArtificialRead(header,"read1",0,1,10);

        ReadPileup pileup = new ReadPileup(null,
                Arrays.asList(read1),
                Arrays.asList(1, 2, 3));
    }


    /**
     * Ensure that splitting read groups still works when dealing with null read groups.
     */
    @Test
    public void testSplitByNullReadGroups() {
        SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(1,1,1000);

        SAMRecord read1 = ArtificialSAMUtils.createArtificialRead(header,"read1",0,1,10);
        SAMRecord read2 = ArtificialSAMUtils.createArtificialRead(header,"read2",0,1,10);
        SAMRecord read3 = ArtificialSAMUtils.createArtificialRead(header,"read3",0,1,10);

        ReadPileup pileup = new ReadPileup(null,
                                                           Arrays.asList(read1, read2, read3),
                                                           Arrays.asList(1, 1, 1));

        ReadPileup nullRgPileup = pileup.getPileupForReadGroup(null);
        List<SAMRecord> nullRgReads = nullRgPileup.getReads();
        Assert.assertEquals(nullRgPileup.size(), 3, "Wrong number of reads in null read group");
        Assert.assertEquals(nullRgReads.get(0), read1, "Read " + read1.getReadName() + " should be in null rg but isn't");
        Assert.assertEquals(nullRgReads.get(1), read2, "Read " + read2.getReadName() + " should be in null rg but isn't");
        Assert.assertEquals(nullRgReads.get(2), read3, "Read " + read3.getReadName() + " should be in null rg but isn't");

        ReadPileup rg1Pileup = pileup.getPileupForReadGroup("rg1");
        Assert.assertNull(rg1Pileup, "Pileup for non-existent read group should return null");
    }

    private static int sampleI = 0;
    private class RBPCountTest {
        final String sample;
        final int nReads, nMapq0, nDeletions;

        private RBPCountTest(int nReads, int nMapq0, int nDeletions) {
            this.sample = "sample" + sampleI++;
            this.nReads = nReads;
            this.nMapq0 = nMapq0;
            this.nDeletions = nDeletions;
        }

        private List<PileupElement> makeReads( final int n, final int mapq, final String op ) {
            final int readLength = 3;

            final List<PileupElement> elts = new LinkedList<>();
            for ( int i = 0; i < n; i++ ) {
                SAMRecord read = ArtificialSAMUtils.createArtificialRead(header, "read", 0, 1, readLength);
                read.setReadBases(Utils.dupBytes((byte) 'A', readLength));
                read.setBaseQualities(Utils.dupBytes((byte) 30, readLength));
                read.setCigarString("1M1" + op + "1M");
                read.setMappingQuality(mapq);
                final int baseOffset = op.equals("M") ? 1 : 0;
                final CigarElement cigarElement = read.getCigar().getCigarElement(1);
                elts.add(new PileupElement(read, baseOffset, cigarElement, 1, 0));
            }

            return elts;
        }

        private ReadPileup makePileup() {
            final List<PileupElement> elts = new LinkedList<>();

            elts.addAll(makeReads(nMapq0, 0, "M"));
            elts.addAll(makeReads(nDeletions, 30, "D"));
            elts.addAll(makeReads(nReads - nMapq0 - nDeletions, 30, "M"));

            return new ReadPileup(loc, elts);
        }

        @Override
        public String toString() {
            return "RBPCountTest{" +
                    "sample='" + sample + '\'' +
                    ", nReads=" + nReads +
                    ", nMapq0=" + nMapq0 +
                    ", nDeletions=" + nDeletions +
                    '}';
        }
    }

    @DataProvider(name = "RBPCountingTest")
    public Object[][] makeRBPCountingTest() {
        final List<Object[]> tests = new LinkedList<>();

        for ( final int nMapq : Arrays.asList(0, 10, 20) ) {
            for ( final int nDeletions : Arrays.asList(0, 10, 20) ) {
                for ( final int nReg : Arrays.asList(0, 10, 20) ) {
                    final int total = nMapq + nDeletions + nReg;
                    if ( total > 0 )
                        tests.add(new Object[]{new RBPCountTest(total, nMapq, nDeletions)});
                }
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "RBPCountingTest")
    public void testRBPCountingTestSinglePileup(RBPCountTest params) {
        testRBPCounts(params.makePileup(), params);
    }

    private void testRBPCounts(final ReadPileup rbp, RBPCountTest expected) {
        for ( int cycles = 0; cycles < 3; cycles++ ) {
            // multiple cycles to make sure caching is working
            Assert.assertEquals(rbp.size(), expected.nReads);
            Assert.assertEquals(rbp.getNumberOfDeletions(), expected.nDeletions);
            Assert.assertEquals(rbp.getNumberOfMappingQualityZeroReads(), expected.nMapq0);
        }
    }

    @Test
    public void testEmptyPileup(){
        final ReadPileup empty = new ReadPileup(loc);
        Assert.assertTrue(empty.isEmpty());
        Assert.assertTrue(empty.size() == 0);
        Assert.assertEquals(empty.getBaseCounts(), new int[4]);
        Assert.assertEquals(empty.getBases(), new byte[0]);
        Assert.assertEquals(empty.getQuals(), new byte[0]);
        Assert.assertEquals(empty.getLocation(), loc);
        Assert.assertEquals(empty.getMappingQuals(), new int[0]);
        Assert.assertTrue(empty.getNegativeStrandPileup().isEmpty());
        Assert.assertTrue(empty.getPositiveStrandPileup().isEmpty());
        Assert.assertTrue(empty.getBaseAndMappingFilteredPileup(0, 0).isEmpty());
        Assert.assertTrue(empty.getBaseFilteredPileup(0).isEmpty());
        Assert.assertEquals(empty.getNumberOfDeletions(), 0);
        Assert.assertEquals(empty.getNumberOfDeletionsAfterThisElement(), 0);
        Assert.assertEquals(empty.getNumberOfInsertionsAfterThisElement(), 0);
        Assert.assertEquals(empty.getNumberOfMappingQualityZeroReads(), 0);
        Assert.assertTrue(empty.getOffsets().isEmpty());
        Assert.assertTrue(empty.getReadGroupIDs().isEmpty());
        Assert.assertTrue(empty.getReads().isEmpty());
        Assert.assertTrue(empty.getSamples().isEmpty());
        Assert.assertTrue(empty.getStartSortedPileup().isEmpty());
        Assert.assertTrue(empty.getPileupForLane("fred") == null);
        Assert.assertTrue(empty.getMappingFilteredPileup(10).isEmpty());
        Assert.assertTrue(empty.toFragments().getOverlappingPairs().isEmpty());
        Assert.assertTrue(empty.toFragments().getSingletonReads().isEmpty());
        Assert.assertEquals("chr1 1 A  ", empty.getPileupString('A'));
    }

    @Test
    public void testSimplePileup(){
        final int readlength = 10;
        final byte[] bases1 = Utils.arrayFromArrayWithLength(new byte[]{'A'}, readlength);
        final byte[] quals1 = Utils.arrayFromArrayWithLength(new byte[]{10}, readlength);
        final String cigar1 = "10M";

        final byte[] bases2 = Utils.arrayFromArrayWithLength(new byte[]{'C'}, readlength);
        final byte[] quals2 = Utils.arrayFromArrayWithLength(new byte[]{20}, readlength);
        final String cigar2 = "5M3I2M";

        final SAMRecord read1 = ArtificialSAMUtils.createArtificialRead(bases1, quals1, cigar1);
        read1.setReadName("read1");
        final SAMRecord read2 = ArtificialSAMUtils.createArtificialRead(bases2, quals2, cigar2);
        read1.setReadName("read2");
        List<SAMRecord> reads = Arrays.asList(read1, read2);
        final ReadPileup pu = new ReadPileup(loc, reads, 1);
        Assert.assertEquals(pu.getBases(), new byte[]{(byte) 'A', (byte) 'C'});
        Assert.assertFalse(pu.isEmpty());
        Assert.assertEquals(pu.size(), 2, "size");
        Assert.assertEquals(pu.getBaseCounts(), new int[]{1,1,0,0});
        Assert.assertEquals(pu.getQuals(), new byte[]{10, 20});
        Assert.assertEquals(pu.getLocation(), loc);
        Assert.assertEquals(pu.getMappingQuals(), new int[]{NO_MAPPING_QUALITY, NO_MAPPING_QUALITY});
        Assert.assertTrue(pu.getNegativeStrandPileup().isEmpty(), "getNegativeStrandPileup");
        Assert.assertEquals(pu.getPositiveStrandPileup().size(), 2, "getPositiveStrandPileup");
        Assert.assertEquals(pu.getBaseAndMappingFilteredPileup(12, 0).size(), 1, "getBaseAndMappingFilteredPileup");
        Assert.assertEquals(pu.getBaseFilteredPileup(12).size(), 1, "getBaseFilteredPileup");
        Assert.assertEquals(pu.getNumberOfDeletions(), 0);
        Assert.assertEquals(pu.getNumberOfDeletionsAfterThisElement(), 0);
        Assert.assertEquals(pu.getNumberOfInsertionsAfterThisElement(), 0);
        Assert.assertEquals(pu.getNumberOfMappingQualityZeroReads(), 2);
        Assert.assertEquals(pu.getOffsets(), Arrays.asList(1, 1), "getOffsets");
        Assert.assertEquals(pu.getReadGroupIDs(), Arrays.asList(new SAMReadGroupRecord[]{null}), "getReadGroups");
        Assert.assertEquals(pu.getReads(), Arrays.asList(read1, read2), "getReads");
        Assert.assertEquals(pu.getSamples(), Arrays.asList(new String[]{null}), "getSamples");
        Assert.assertEquals(pu.getStartSortedPileup().size(), 2, "getStartSortedPileup");
        Assert.assertTrue(pu.getPileupForLane("fred") == null);
        Assert.assertTrue(pu.getMappingFilteredPileup(10).isEmpty());
        Assert.assertTrue(pu.toFragments().getOverlappingPairs().isEmpty());
        Assert.assertEquals(pu.toFragments().getSingletonReads().stream().map(pe -> pe.getRead()).collect(Collectors.toList()), Arrays.asList(read1, read2), "getSingletonReads");
        Assert.assertEquals("chr1 1 A AC +5", pu.getPileupString('A'));
    }

    @Test
    public void testSimplePileupWithOffset(){
        final int readlength = 10;
        final byte[] bases1 = Utils.arrayFromArrayWithLength(new byte[]{'A'}, readlength);
        final byte[] quals1 = Utils.arrayFromArrayWithLength(new byte[]{10}, readlength);
        final String cigar1 = "10M";

        final byte[] bases2 = Utils.arrayFromArrayWithLength(new byte[]{'C'}, readlength);
        final byte[] quals2 = Utils.arrayFromArrayWithLength(new byte[]{20}, readlength);
        final String cigar2 = "10M";

        final SAMRecord read1 = ArtificialSAMUtils.createArtificialRead(bases1, quals1, cigar1);
        read1.setReadName("read1");
        final SAMRecord read2 = ArtificialSAMUtils.createArtificialRead(bases2, quals2, cigar2);
        read1.setReadName("read2");
        List<SAMRecord> reads = Arrays.asList(read1, read2);
        final int off = 6;
        final ReadPileup pu = new ReadPileup(loc, reads, off);
        Assert.assertEquals(pu.getBases(), new byte[]{(byte) 'A', (byte) 'C'});
        Assert.assertFalse(pu.isEmpty());
        Assert.assertEquals(pu.size(), 2, "size");
        Assert.assertEquals(pu.getBaseCounts(), new int[]{1,1,0,0});
        Assert.assertEquals(pu.getQuals(), new byte[]{10, 20});
        Assert.assertEquals(pu.getLocation(), loc);
        Assert.assertEquals(pu.getMappingQuals(), new int[]{NO_MAPPING_QUALITY, NO_MAPPING_QUALITY});
        Assert.assertTrue(pu.getNegativeStrandPileup().isEmpty(), "getNegativeStrandPileup");
        Assert.assertEquals(pu.getPositiveStrandPileup().size(), 2, "getPositiveStrandPileup");
        Assert.assertEquals(pu.getBaseAndMappingFilteredPileup(12, 0).size(), 1, "getBaseAndMappingFilteredPileup");
        Assert.assertEquals(pu.getBaseFilteredPileup(12).size(), 1, "getBaseFilteredPileup");
        Assert.assertEquals(pu.getNumberOfDeletions(), 0);
        Assert.assertEquals(pu.getNumberOfDeletionsAfterThisElement(), 0);
        Assert.assertEquals(pu.getNumberOfInsertionsAfterThisElement(), 0);
        Assert.assertEquals(pu.getNumberOfMappingQualityZeroReads(), 2);
        Assert.assertEquals(pu.getOffsets(), Arrays.asList(off, off), "getOffsets");
        Assert.assertEquals(pu.getReadGroupIDs(), Arrays.asList(new SAMReadGroupRecord[]{null}), "getReadGroups");
        Assert.assertEquals(pu.getReads(), Arrays.asList(read1, read2), "getReads");
        Assert.assertEquals(pu.getSamples(), Arrays.asList(new String[]{null}), "getSamples");
        Assert.assertEquals(pu.getStartSortedPileup().size(), 2, "getStartSortedPileup");
        Assert.assertTrue(pu.getPileupForLane("fred") == null);
        Assert.assertTrue(pu.getMappingFilteredPileup(10).isEmpty());
        Assert.assertTrue(pu.toFragments().getOverlappingPairs().isEmpty());
        Assert.assertEquals(pu.toFragments().getSingletonReads().stream().map(pe -> pe.getRead()).collect(Collectors.toList()), Arrays.asList(read1, read2), "getSingletonReads");
        Assert.assertEquals("chr1 1 A AC +5", pu.getPileupString('A'));
    }
}