package org.broadinstitute.hellbender.utils.pileup;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

import static org.broadinstitute.hellbender.utils.read.ReadConstants.NO_MAPPING_QUALITY;

/**
 * Test routines for read-backed pileup.
 */
public final class ReadPileupUnitTest extends GATKBaseTest {
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

    @Test
    public void testGetPileupForSample() {
        // read groups and header
        final SAMReadGroupRecord[] readGroups = new SAMReadGroupRecord[5];
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader(1, 1, 1000);
        int s = 0;
        // readGroups[4] is left null intentionally
        for(final String sample: Arrays.asList("sample1", "sample1", "sample2", null)) {
            readGroups[s] = new SAMReadGroupRecord("rg"+s);
            readGroups[s].setSample(sample);
            header.addReadGroup(readGroups[s++]);
        }

        // reads
        final int rgCoverage = 4;
        final List<GATKRead> reads = new ArrayList<>(rgCoverage*readGroups.length);
        for(int i = 0; i < rgCoverage; i++) {
            for(final SAMReadGroupRecord rg: readGroups) {
                final GATKRead r = ArtificialReadUtils.createArtificialRead(header, (rg == null) ? "null" : rg.getReadGroupId() + "_" + rg.getSample() + "_" + i, 0, 1, 10);
                if(rg != null) {
                    r.setReadGroup(rg.getId());
                }
                reads.add(r);
            }
        }

        // pileup
        final ReadPileup pileup = new ReadPileup(loc, reads, 1);
        // sample1
        final ReadPileup pileupSample1 = pileup.getPileupForSample("sample1", header);
        Assert.assertEquals(pileupSample1.size(), rgCoverage*2, "Wrong number of elements for sample1");
        Assert.assertTrue( pileupSample1.getReadGroupIDs().contains("rg0"), "Pileup for sample1 should contain rg0");
        Assert.assertTrue( pileupSample1.getReadGroupIDs().contains("rg1"), "Pileup for sample1 should contain rg1");
        Assert.assertFalse(pileupSample1.getReadGroupIDs().contains("rg2"), "Pileup for sample1 shouldn't contain rg2");
        Assert.assertFalse(pileupSample1.getReadGroupIDs().contains("rg3"), "Pileup for sample1 shouldn't contain rg3");
        Assert.assertFalse(pileupSample1.getReadGroupIDs().contains(null),  "Pileup for sample1 shouldn't contain null read group");

        // sample2
        final ReadPileup pileupSample2 = pileup.getPileupForSample("sample2", header);
        Assert.assertEquals(pileupSample2.size(), rgCoverage, "Wrong number of elements for sample2");
        Assert.assertFalse(pileupSample2.getReadGroupIDs().contains("rg0"), "Pileup for sample2 shouldn't contain rg0");
        Assert.assertFalse(pileupSample2.getReadGroupIDs().contains("rg1"), "Pileup for sample2 shouldn't contain rg1");
        Assert.assertTrue( pileupSample2.getReadGroupIDs().contains("rg2"), "Pileup for sample2 should contain rg2");
        Assert.assertFalse(pileupSample2.getReadGroupIDs().contains("rg3"), "Pileup for sample2 shouldn't contain rg3");
        Assert.assertFalse(pileupSample2.getReadGroupIDs().contains(null),  "Pileup for sample2 shouldn't contain null read group");

        // null sample
        final ReadPileup pileupNull = pileup.getPileupForSample(null, header);
        Assert.assertEquals(pileupNull.size(), rgCoverage*2, "Wrong number of elements for null sample");
        Assert.assertFalse(pileupNull.getReadGroupIDs().contains("rg0"), "Pileup for null sample shouldn't contain rg0");
        Assert.assertFalse(pileupNull.getReadGroupIDs().contains("rg1"), "Pileup for null sample shouldn't contain rg1");
        Assert.assertFalse(pileupNull.getReadGroupIDs().contains("rg2"), "Pileup for null sample shouldn't contain rg2");
        Assert.assertTrue( pileupNull.getReadGroupIDs().contains("rg3"), "Pileup for null sample should contain rg3");
        Assert.assertTrue( pileupNull.getReadGroupIDs().contains(null),  "Pileup for null sample should contain null read group");

    }

    @Test(expectedExceptions = UserException.ReadMissingReadGroup.class)
    public void testSplitBySampleMissingReadGroupException() throws Exception {
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader(1, 1, 1000);
        final GATKRead read = ArtificialReadUtils.createArtificialRead(header, "nullRG", 0, 1, 10);
        final ReadPileup pileup = new ReadPileup(loc, Collections.singletonList(read), 1);
        final Map<String, ReadPileup> error = pileup.splitBySample(header, null);
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
        Assert.assertEquals(pu.getPileupString('N'), String.format("%s %s N %s%s %s%s",
                loc.getContig(), loc.getStart(), // the position
                (char) read1.getBase(0), (char) read2.getBase(0),  // the bases
                SAMUtils.phredToFastq(read1.getBaseQuality(0)), // base quality in FASTQ format
                SAMUtils.phredToFastq(read2.getBaseQuality(0)))); // base quality in FASTQ format

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

        final PileupElement firstElem = pu.getElementForRead(read1);
        Assert.assertNull(firstElem.getAdjacentOperator(PileupElement.Direction.NEXT));
        Assert.assertNull(firstElem.getAdjacentOperator(PileupElement.Direction.PREV));
        final PileupElement secondElem = pu.getElementForRead(read2);
        Assert.assertEquals(secondElem.getAdjacentOperator(PileupElement.Direction.NEXT), CigarOperator.I);
        Assert.assertNull(secondElem.getAdjacentOperator(PileupElement.Direction.PREV));

        Assert.assertNotNull(pu.toString()); //checking non-blowup. We're not making any claims about the format of toString
        Assert.assertEquals(pu.getPileupString('N'), String.format("%s %s N %s%s %s%s",
                loc.getContig(), loc.getStart(), // the position
                (char) read1.getBase(0), (char) read2.getBase(0),  // the bases
                SAMUtils.phredToFastq(read1.getBaseQuality(0)), // base quality in FASTQ format
                SAMUtils.phredToFastq(read2.getBaseQuality(0)))); // base quality in FASTQ format

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
        Assert.assertEquals(pu.getPileupString('N'), String.format("%s %s N %s%s %s%s",
                loc.getContig(), loc.getStart(), // the position
                (char) read1.getBase(off), (char) read2.getBase(off),  // the bases
                SAMUtils.phredToFastq(read1.getBaseQuality(off)), // base quality in FASTQ format
                SAMUtils.phredToFastq(read2.getBaseQuality(off)))); // base quality in FASTQ format

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

    @DataProvider(name = "FixPairOverlappingQualitiesTest")
    public Object[][] overlappingElementsToFix() throws Exception {
        final byte highQual = 60;
        final byte lowQual = 10;
        final byte qualitySum = 70;
        final byte reducedQuality = 48; // 80% of 60
        final byte zeroQuality = 0;
        final GATKRead read1 = ArtificialReadUtils.createArtificialRead("6M");
        final GATKRead read2 = ArtificialReadUtils.createArtificialRead("6M");
        // set the paired and mate state
        read1.setIsPaired(true);
        read1.setMatePosition(read2);
        // set bases and qualities
        read1.setBases(new byte[] {'A', 'A', 'A', 'A', 'T', 'T'});
        read1.setBaseQualities(new byte[] {highQual, highQual, highQual, lowQual,
                highQual, SAMUtils.MAX_PHRED_SCORE});
        // set the paired and mate state
        read2.setIsPaired(true);
        read2.setMatePosition(read1);
        // set bases and qualities
        read2.setBases(new byte[] {'A', 'T', 'T', 'T', 'T', 'T'});
        read2.setBaseQualities(new byte[] {lowQual, highQual, lowQual, highQual,
                highQual, SAMUtils.MAX_PHRED_SCORE});
        return new Object[][] {
                // Same base, first element with higher quality
                {PileupElement.createPileupForReadAndOffset(read1, 0),
                        PileupElement.createPileupForReadAndOffset(read2, 0),
                        qualitySum, zeroQuality},
                // Different base, both with same quality
                {PileupElement.createPileupForReadAndOffset(read1, 1),
                        PileupElement.createPileupForReadAndOffset(read2, 1),
                        reducedQuality, zeroQuality},
                // Different base, first with higher quality
                {PileupElement.createPileupForReadAndOffset(read1, 2),
                        PileupElement.createPileupForReadAndOffset(read2, 2),
                        reducedQuality, zeroQuality},
                // Different base, second with higher quality
                {PileupElement.createPileupForReadAndOffset(read1, 3),
                        PileupElement.createPileupForReadAndOffset(read2, 3),
                        zeroQuality, reducedQuality},
                // Same base, high quality for both (simple cap)
                {PileupElement.createPileupForReadAndOffset(read1, 4),
                        PileupElement.createPileupForReadAndOffset(read2, 4),
                        QualityUtils.MAX_SAM_QUAL_SCORE, zeroQuality},
                // Same base, maximum qualities for both (simple cap)
                {PileupElement.createPileupForReadAndOffset(read1, 5),
                        PileupElement.createPileupForReadAndOffset(read2, 5),
                        QualityUtils.MAX_SAM_QUAL_SCORE, zeroQuality},
        };
    }

    @Test(dataProvider = "FixPairOverlappingQualitiesTest")
    public void testFixPairOverlappingQualities(final PileupElement first,
                                                final PileupElement second, final byte expectedQualFirst,
                                                final byte expectedQualSecond) throws Exception {
        ReadPileup.fixPairOverlappingQualities(first, second);
        Assert.assertEquals(first.getQual(), expectedQualFirst);
        Assert.assertEquals(second.getQual(), expectedQualSecond);
    }


    @Test(dataProvider = "FixPairOverlappingQualitiesTest")
    public void testFixPairOverlappingQualitiesFromReadPileup(final PileupElement first,
            final PileupElement second, final byte expectedQualFirst,
            final byte expectedQualSecond) throws Exception {
        final List<PileupElement> elements = Arrays.asList(first, second);
        final ReadPileup pileup = new ReadPileup(loc, elements);
        pileup.fixOverlaps();
        Assert.assertEquals(elements.get(0).getQual(), expectedQualFirst);
        Assert.assertEquals(elements.get(1).getQual(), expectedQualSecond);
    }

    @Test
    public void testFixPairOverlappingQualitiesCap() {
        final PileupElement element1 = PileupElement
                .createPileupForReadAndOffset(ArtificialReadUtils.createArtificialRead("1M"), 0);
        element1.getRead().setName("read1");
        element1.getRead().setBases(new byte[] {'A'});
        final PileupElement element2 = PileupElement
                .createPileupForReadAndOffset(ArtificialReadUtils.createArtificialRead("1M"), 0);
        element2.getRead().setName("read2");
        element2.getRead().setBases(new byte[] {'A'});
        final PileupElement element3 = PileupElement
                .createPileupForReadAndOffset(ArtificialReadUtils.createArtificialRead("1M"), 0);
        element3.getRead().setName("read3");
        element2.getRead().setBases(new byte[] {'A'});
        element3.getRead().setBaseQualities(new byte[] {Byte.MAX_VALUE});
        // Check all possible combinations that goes beyond the maximum value
        for (byte i = 0; i < Byte.MAX_VALUE; i++) {
            final byte[] iArray = new byte[] {i};
            final byte j = (byte) (Byte.MAX_VALUE - i);
            element1.getRead().setBaseQualities(iArray);
            element2.getRead().setBaseQualities(new byte[] {j});
            logger.debug("Test: fixing ({}) and ({})", element1, element2);
            ReadPileup.fixPairOverlappingQualities(element1, element2);
            Assert.assertEquals(element1.getQual(), QualityUtils.MAX_SAM_QUAL_SCORE);
            Assert.assertEquals(element2.getQual(), 0);
            element1.getRead().setBaseQualities(iArray);
            logger.debug("Test: fixing ({}) and ({})", element3, element1);
            ReadPileup.fixPairOverlappingQualities(element3, element1);
            Assert.assertEquals(element3.getQual(), QualityUtils.MAX_SAM_QUAL_SCORE);
            Assert.assertEquals(element1.getQual(), 0);
        }
        // Finally, check what happens if both are the maximum value
        element1.getRead().setBaseQualities(new byte[] {Byte.MAX_VALUE});
        logger.debug("Test: fixing ({}) and ({})", element3, element1);
        ReadPileup.fixPairOverlappingQualities(element3, element1);
        Assert.assertEquals(element3.getQual(), QualityUtils.MAX_SAM_QUAL_SCORE);
        Assert.assertEquals(element1.getQual(), 0);
    }

    @Test
    public void testFixOverlappingQualitiesDeletedElements() throws Exception {
        // Create two reads with deletion
        final GATKRead read = ArtificialReadUtils.createArtificialRead("1M1D1M");
        read.setBases(new byte[] {'A', 'A'});
        read.setBaseQualities(new byte[] {60, 60});
        final GATKRead read2 = ArtificialReadUtils.createArtificialRead("1M1D1M");
        read.setBases(new byte[] {'A', 'A'});
        read2.setBaseQualities(new byte[] {40, 40});
        // Creates PileupElement with a copy of the reads to keep the original read unmodified, because the read from
        // the PileupElement is a reference to the Object passed to the constructor. This is required
        // because fixing the overlapping qualities change the read stored into PileupElement and we have to check that
        // it does not change anything from the read.
        final GATKRead copy = read.deepCopy();
        final GATKRead copy2 = read2.deepCopy();
        final PileupElement deleted = new PileupElement(copy, 1, copy.getCigarElement(1), 1, 0);
        final PileupElement normal = new PileupElement(copy2, 0, copy.getCigarElement(0), 0, 0);
        ReadPileup.fixPairOverlappingQualities(deleted, normal);
        // Copy reads are a reference to the read into the PileupElement, and thus should not be modified.
        // Thus, we use them to test if it is not modified. This could be checked uncommenting the next two lines:
        // Assert.assertSame(deleted.getRead(), copy);
        // Assert.assertSame(normal.getRead(), copy2);
        Assert.assertEquals(copy, read);
        Assert.assertEquals(copy2, read2);
        ReadPileup.fixPairOverlappingQualities(normal, deleted);
        Assert.assertEquals(copy, read);
        Assert.assertEquals(copy2, read2);
    }

    @Test
    public void testSortedIterator() throws Exception {
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader();
        final GATKRead read1 = ArtificialReadUtils.createArtificialRead(header, "10M");
        final GATKRead read2 = ArtificialReadUtils.createArtificialRead(header, "10M");
        read1.setPosition(new SimpleInterval("22", 200, 210));
        read2.setPosition(new SimpleInterval("22", 208, 218));
        read1.setMatePosition(read2);
        read2.setMatePosition(read1);
        final Locatable loc = new SimpleInterval("22", 208, 208);
        final Map<String, ReadPileup> stratified = new LinkedHashMap<>();
        stratified.put("sample1", new ReadPileup(loc, Arrays.asList(read2), 0));
        stratified.put("sample2", new ReadPileup(loc, Arrays.asList(read1), 9));
        final ReadPileup combined = new ReadPileup(loc, stratified);
        final Iterator<PileupElement> sortedIterator = combined.sortedIterator();
        Assert.assertSame(sortedIterator.next().getRead(), read1);
        Assert.assertSame(sortedIterator.next().getRead(), read2);
    }

    @Test
    public void testOverlappingFragmentFilter() throws Exception {
        final String rgID = "MY.ID";
        final SAMReadGroupRecord rg = new SAMReadGroupRecord(rgID);
        rg.setSample("MySample");
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeaderWithReadGroup(rg);
        final int readlength = 10;
        final String cigar = "10M";
        final byte[] lowQuals = Utils.repeatBytes((byte)10, readlength);
        final byte[] highQuals = Utils.repeatBytes((byte)20, readlength);
        final byte[] basesAllA = Utils.repeatChars('A', readlength);
        final byte[] basesAllC = Utils.repeatChars('C', readlength);
        final byte[] basesAllT = Utils.repeatChars('T', readlength);

        GATKRead read1BasesDisagree = ArtificialReadUtils.createArtificialRead(header, basesAllA, lowQuals, cigar);
        GATKRead read2BasesDisagree = ArtificialReadUtils.createArtificialRead(header, basesAllC, highQuals, cigar);
        GATKRead read3BasesDisagree = ArtificialReadUtils.createArtificialRead(header, basesAllT, lowQuals, cigar);
        read1BasesDisagree.setName("BasesDisagree");
        read2BasesDisagree.setName("BasesDisagree");
        read3BasesDisagree.setName("BasesDisagree");
        read1BasesDisagree.setReadGroup(rgID);
        read2BasesDisagree.setReadGroup(rgID);
        read3BasesDisagree.setReadGroup(rgID);

        GATKRead read1BasesAgree = ArtificialReadUtils.createArtificialRead(header, basesAllA, lowQuals, cigar);
        GATKRead read2BasesAgree = ArtificialReadUtils.createArtificialRead(header, basesAllA, highQuals, cigar);
        read1BasesAgree.setName("BasesAgree");
        read2BasesAgree.setName("BasesAgree");
        read1BasesAgree.setReadGroup(rgID);
        read2BasesAgree.setReadGroup(rgID);

        final List<GATKRead> reads = Arrays.asList(read1BasesDisagree, read2BasesDisagree, read1BasesAgree, read2BasesAgree);
        final int off = 0;
        final ReadPileup pu = new ReadPileup(loc, reads, off);

        final ReadPileup filteredPileupDiscardDiscordant = pu.getOverlappingFragmentFilteredPileup(header);
        final List<GATKRead> filteredReadsDiscardDiscordant = filteredPileupDiscardDiscordant.getReads();
        Assert.assertFalse(filteredReadsDiscardDiscordant.contains(read1BasesDisagree), "Reads with disagreeing bases were kept.");
        Assert.assertFalse(filteredReadsDiscardDiscordant.contains(read2BasesDisagree), "Reads with disagreeing bases were kept.");
        Assert.assertFalse(filteredReadsDiscardDiscordant.contains(read3BasesDisagree), "Reads with disagreeing bases were kept.");
        Assert.assertFalse(filteredReadsDiscardDiscordant.contains(read1BasesAgree), "The lower quality base was kept.");
        Assert.assertTrue(filteredReadsDiscardDiscordant.contains(read2BasesAgree), "The higher quality base is missing.");

        final ReadPileup filteredPileupKeepDiscordant = pu.getOverlappingFragmentFilteredPileup(false, ReadPileup.baseQualTieBreaker, header);
        final List<GATKRead> filteredReadsKeepDiscordant = filteredPileupKeepDiscordant.getReads();
        Assert.assertFalse(filteredReadsKeepDiscordant.contains(read1BasesDisagree), "Low quality base was kept.");
        Assert.assertTrue(filteredReadsKeepDiscordant.contains(read2BasesDisagree), "High quality base is missing.");
        Assert.assertFalse(filteredReadsKeepDiscordant.contains(read3BasesDisagree), "Low quality base was kept.");
        Assert.assertFalse(filteredReadsKeepDiscordant.contains(read1BasesAgree), "The lower quality base was kept.");
        Assert.assertTrue(filteredReadsKeepDiscordant.contains(read2BasesAgree), "The higher quality base is missing.");
    }
}