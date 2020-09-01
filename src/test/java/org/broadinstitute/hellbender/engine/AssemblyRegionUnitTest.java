package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.Assert;
import org.testng.annotations.AfterClass;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;


public final class AssemblyRegionUnitTest extends GATKBaseTest {
    private ReferenceSequenceFile seq;
    private String contig;
    private int contigLength;
    private SAMFileHeader header;

    @BeforeClass
    public void init() {
        // sequence
        seq = new CachingIndexedFastaSequenceFile(IOUtils.getPath(hg19MiniReference));
        contig = "1";
        contigLength = seq.getSequence(contig).length();
        header = ArtificialReadUtils.createArtificialSamHeader(seq.getSequenceDictionary());
    }

    @AfterClass
    public void close() throws IOException {
        seq.close();
    }

    @Test
    public void testConstructor(){
        final SimpleInterval loc = new SimpleInterval("1", 10, 20);
        final AssemblyRegion ar = new AssemblyRegion(loc, 2, header);
        Assert.assertEquals(ar.isActive(), true);
        Assert.assertEquals(ar.getSpan(), loc);
        Assert.assertEquals(ar.getHeader(), header);
    }

    @DataProvider(name = "ActionRegionCreationTest")
    public Object[][] makePollingData() {
        final List<Object[]> tests = new ArrayList<>();
        for ( final int start : Arrays.asList(1, 10, 100, contigLength - 10, contigLength - 1) ) {
            for ( final int size : Arrays.asList(1, 10, 100, 1000) ) {
                for ( final int ext : Arrays.asList(0, 1, 10, 100) ) {
                    for ( final boolean isActive : Arrays.asList(true, false) ) {
                        final SimpleInterval loc = IntervalUtils.trimIntervalToContig(contig, start, start + size - 1, contigLength);
                        tests.add(new Object[]{loc, isActive, ext});
                    }
                }
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "ActionRegionCreationTest")
    public void testCreatingAssemblyRegions(final SimpleInterval loc, final boolean isActive, final int extension) {
        final AssemblyRegion region = new AssemblyRegion(loc, isActive, extension, header);
        Assert.assertFalse(region.isFinalized());
        Assert.assertEquals(region.getSpan(), loc);
        Assert.assertEquals(region.getPaddedSpan().getStart(), Math.max(loc.getStart() - extension, 1));
        Assert.assertEquals(region.getPaddedSpan().getEnd(), Math.min(loc.getEnd() + extension, contigLength));
        Assert.assertEquals(region.isActive(), isActive);
        Assert.assertEquals(region.getReads(), Collections.emptyList());
        Assert.assertEquals(region.size(), 0);
        Assert.assertNotNull(region.toString());

        assertGoodReferenceGetter(region.getAssemblyRegionReference(seq), region.getPaddedSpan(), 0);
        assertGoodReferenceGetter(region.getAssemblyRegionReference(seq, 0), region.getPaddedSpan(), 0);
        assertGoodReferenceGetter(region.getAssemblyRegionReference(seq, 10), region.getPaddedSpan(), 10);

        region.setFinalized(false);
        Assert.assertFalse(region.isFinalized());
        region.setFinalized(true);
        Assert.assertTrue(region.isFinalized());
        region.setFinalized(false);
        Assert.assertFalse(region.isFinalized());
    }

    private void assertGoodReferenceGetter(final byte[] actualBytes, final SimpleInterval span, final int padding) {
        final int expectedStart = Math.max(span.getStart() - padding, 1);
        final int expectedStop = Math.min(span.getEnd() + padding, contigLength);
        final byte[] expectedBytes = seq.getSubsequenceAt(span.getContig(), expectedStart, expectedStop).getBases();
        Assert.assertEquals(actualBytes, expectedBytes);
    }

    @DataProvider(name = "AssemblyRegionReads")
    public Object[][] makeAssemblyRegionReads() {
        final List<Object[]> tests = new ArrayList<>();
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader(seq.getSequenceDictionary());
        for ( final int start : Arrays.asList(1, 10, 100, contigLength - 10, contigLength - 1) ) {
            for ( final int readStartOffset : Arrays.asList(-100, -10, 0, 10, 100) ) {
                for ( final int readSize : Arrays.asList(10, 100, 1000) ) {
                    final SimpleInterval loc = IntervalUtils.trimIntervalToContig(contig, start, start + 10, header.getSequence(contig).getSequenceLength());

                    final int readStart = Math.max(start + readStartOffset, 1);
                    final int readStop = Math.min(readStart + readSize, contigLength);
                    final int readLength = readStop - readStart + 1;
                    if ( readLength > 0 ) {
                        final GATKRead read = ArtificialReadUtils.createArtificialRead(header, "read", 0, readStart, readLength);
                        final SimpleInterval readLoc = new SimpleInterval(read);
                        if ( readLoc.overlaps(loc) ) {
                            tests.add(new Object[]{loc, read});
                        }
                    }
                }
            }
        }

        return tests.subList(2,3).toArray(new Object[][]{});       //HACK!
    }

    @Test(expectedExceptions = UserException.MissingContigInSequenceDictionary.class)
    //Testing a case where an NPE was thrown when somehow the reference sequence dictionary didn't contain a contig in the
    //assembly region.
    public void testMismatchedReferenceAndRegion(){
        final SAMSequenceDictionary sequenceDictionary = new SAMSequenceDictionary(new ArrayList<>(header.getSequenceDictionary().getSequences()));
        final String contigNotInReference = "chrNotInReference";
        sequenceDictionary.addSequence(new SAMSequenceRecord(contigNotInReference, 1000));
        //it's not possible to create an assembly region on a contig that isn't in the given SamHeader
        final SAMFileHeader modifiedHeader = new SAMFileHeader(sequenceDictionary);
        final AssemblyRegion assemblyRegion = new AssemblyRegion(new SimpleInterval(contigNotInReference,1, 200), 10, modifiedHeader);
        assemblyRegion.getAssemblyRegionReference(hg19ReferenceReader, 10);
    }

    @Test(dataProvider = "AssemblyRegionReads")
    public void testAssemblyRegionReads(final SimpleInterval loc, final GATKRead read) throws Exception {
        final AssemblyRegion region = new AssemblyRegion(loc, true, 0, header);
        Assert.assertEquals(region.getReads(), Collections.emptyList());
        Assert.assertEquals(region.size(), 0);
        Assert.assertEquals(region.getPaddedSpan(), loc);

        region.add(read);
        Assert.assertEquals(region.getReads(), Collections.singletonList(read));
        Assert.assertEquals(region.size(), 1);
        Assert.assertEquals(region.getPaddedSpan(), loc);

        region.clearReads();
        Assert.assertEquals(region.getReads(), Collections.emptyList());
        Assert.assertEquals(region.size(), 0);
        Assert.assertEquals(region.getPaddedSpan(), loc);

        region.addAll(Collections.singleton(read));
        Assert.assertEquals(region.getReads(), Collections.singletonList(read));
        Assert.assertEquals(region.size(), 1);
        Assert.assertEquals(region.getPaddedSpan(), loc);

        region.removeAll(Collections.<GATKRead>emptySet());
        Assert.assertEquals(region.getReads(), Collections.singletonList(read));
        Assert.assertEquals(region.size(), 1);
        Assert.assertEquals(region.getPaddedSpan(), loc);

        region.removeAll(Collections.singleton(read));
        Assert.assertEquals(region.getReads(), Collections.emptyList());
        Assert.assertEquals(region.size(), 0);
        Assert.assertEquals(region.getPaddedSpan(), loc);

        final GATKRead read2 = read.copy();
        read2.setName(read.getName() + ".clone");

        for ( final GATKRead readToKeep : Arrays.asList(read, read2)) {
            region.addAll(Arrays.asList(read, read2));
            final GATKRead readToDiscard = readToKeep == read ? read2 : read;
            region.removeAll(Collections.singleton(readToDiscard));
            Assert.assertEquals(region.getReads(), Arrays.asList(readToKeep));
            Assert.assertEquals(region.size(), 1);
            Assert.assertEquals(region.getPaddedSpan(), loc);
        }
    }

    // -----------------------------------------------------------------------------------------------
    //
    // Make sure bad inputs are properly detected
    //
    // -----------------------------------------------------------------------------------------------

    @DataProvider(name = "BadReadsTest")
    public Object[][] makeBadReadsTest() {
        List<Object[]> tests = new ArrayList<>();
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader(seq.getSequenceDictionary());
        tests.add(new Object[]{
                header,
                ArtificialReadUtils.createArtificialRead(header, "read1", 0, 10, 10),
                ArtificialReadUtils.createArtificialRead(header, "read2", 0, 9, 10)});
        tests.add(new Object[]{
                header,
                ArtificialReadUtils.createArtificialRead(header, "read1", 0, 10, 10),
                ArtificialReadUtils.createArtificialRead(header, "read2", 1, 9, 10)});
        tests.add(new Object[]{
                header,
                ArtificialReadUtils.createArtificialRead(header, "read1", 1, 10, 10),
                ArtificialReadUtils.createArtificialRead(header, "read2", 0, 9, 10)});
        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "BadReadsTest", expectedExceptions = IllegalArgumentException.class)
    public void testBadReads(final SAMFileHeader header, final GATKRead read1, final GATKRead read2) {
        final SimpleInterval loc = new SimpleInterval(read1);
        final AssemblyRegion region = new AssemblyRegion(loc, true, 0, header);
        region.add(read1);
        region.add(read2);
    }

    // -----------------------------------------------------------------------------------------------
    //
    // Make sure we can properly cut up an assembly region based on engine intervals
    //
    // -----------------------------------------------------------------------------------------------

    @DataProvider(name = "TrimAssemblyRegionData")
    public Object[][] makeTrimAssemblyRegionData() {
        final List<Object[]> tests = new ArrayList<>();

        // fully enclosed within active region
        tests.add(new Object[]{
                new SimpleInterval("1", 10, 20), 10,
                new SimpleInterval("1", 15, 16),
                new SimpleInterval("1", 15, 16), 0});

        tests.add(new Object[]{
                new SimpleInterval("1", 10, 20), 10,
                new SimpleInterval("1", 10, 15),
                new SimpleInterval("1", 10, 15), 0});

        tests.add(new Object[]{
                new SimpleInterval("1", 10, 20), 10,
                new SimpleInterval("1", 15, 20),
                new SimpleInterval("1", 15, 20), 0});

        // needs extra padding on the right
        tests.add(new Object[]{
                new SimpleInterval("1", 10, 20), 10,
                new SimpleInterval("1", 15, 25),
                new SimpleInterval("1", 15, 20), 5});

        // needs extra padding on the left
        tests.add(new Object[]{
                new SimpleInterval("1", 10, 20), 10,
                new SimpleInterval("1", 5, 15),
                new SimpleInterval("1", 10, 15), 5});

        // needs extra padding on both
        tests.add(new Object[]{
                new SimpleInterval("1", 10, 20), 10,
                new SimpleInterval("1", 7, 21),
                new SimpleInterval("1", 10, 20), 3});
        tests.add(new Object[]{
                new SimpleInterval("1", 10, 20), 10,
                new SimpleInterval("1", 9, 23),
                new SimpleInterval("1", 10, 20), 3});

        // desired span captures everything, so we're returning everything.  Tests that extension is set correctly
        tests.add(new Object[]{
                new SimpleInterval("1", 10, 20), 10,
                new SimpleInterval("1", 1, 50),
                new SimpleInterval("1", 10, 20), 10});

        // At the start of the chromosome, potentially a bit weird
        tests.add(new Object[]{
                new SimpleInterval("1", 1, 10), 10,
                new SimpleInterval("1", 1, 50),
                new SimpleInterval("1", 1, 10), 10});

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "TrimAssemblyRegionData")
    public void testTrimAssemblyRegion(final SimpleInterval regionLoc, final int extension, final SimpleInterval desiredSpan, final SimpleInterval expectedAssemblyRegion, final int expectedExtension) {
        final AssemblyRegion region = new AssemblyRegion(regionLoc, true, extension, header);
        final AssemblyRegion trimmed = region.trim(desiredSpan, desiredSpan);
        Assert.assertEquals(trimmed.getSpan(), expectedAssemblyRegion, "Incorrect region");
    }
}
