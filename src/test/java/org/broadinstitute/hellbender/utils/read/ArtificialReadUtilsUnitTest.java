package org.broadinstitute.hellbender.utils.read;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.TextCigarCodec;
import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.utils.pileup.PileupElement;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Collection;
import java.util.Iterator;

import static org.testng.Assert.*;

public final class ArtificialReadUtilsUnitTest extends GATKBaseTest {

    @Test
    public void testReadGroup() {
        GATKRead read = ArtificialReadUtils.createArtificialRead(TextCigarCodec.decode("10M"));
        assertNull(read.getReadGroup());
    }

    @Test
    public void testCreateIdenticalArtificialReads(){
        int length= 7;
        int alignmentStart= 19;
        int refIndex= 0;
        SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader();
        final Collection<GATKRead> reads = ArtificialReadUtils.createIdenticalArtificialReads(10, header, "name", refIndex, alignmentStart, length);
        Assert.assertEquals(reads.size(), 10);
        for (final GATKRead read: reads){
            Assert.assertEquals(read.getLength(), length);
            Assert.assertEquals(read.getStart(), alignmentStart);
            Assert.assertEquals(ReadUtils.getReferenceIndex(read, header), refIndex);
        }
    }

    @Test
    public void basicReadIteratorTest() {
        Iterator<GATKRead> iter = ArtificialReadUtils.mappedReadIterator(1, 100, 100);
        int count = 0;
        while (iter.hasNext()) {
            GATKRead rec = iter.next();
            count++;
        }
        assertEquals(count, 100 * 100);
    }

    @Test
    public void tenPerChromosome() {
        ArtificialReadQueryIterator iter = ArtificialReadUtils.mappedReadIterator(1, 100, 10);
        int count = 0;
        while (iter.hasNext()) {
            GATKRead rec = iter.next();

            assertEquals(Math.round(count / 10), ReadUtils.getReferenceIndex(rec, iter.getHeader()));
            count++;
        }
        assertEquals(count, 100 * 10);
    }

    @Test
    public void onePerChromosome() {
        ArtificialReadQueryIterator iter = ArtificialReadUtils.mappedReadIterator(1, 100, 1);
        int count = 0;
        while (iter.hasNext()) {
            GATKRead rec = iter.next();

            assertEquals(count, ReadUtils.getReferenceIndex(rec, iter.getHeader()));
            count++;
        }
        assertEquals(count, 100);
    }

    @Test
    public void basicUnmappedIteratorTest() {
        Iterator<GATKRead> iter = ArtificialReadUtils.mappedAndUnmappedReadIterator(1, 100, 100, 1000);
        int count = 0;
        for (int x = 0; x < (100* 100); x++ ) {
            if (!iter.hasNext()) {
                fail ("we didn't get the expected number of reads");
            }
            GATKRead rec = iter.next();
            assertFalse(rec.isUnmapped());
            count++;
        }
        assertEquals(100 * 100, count);

        // now we should have 1000 unmapped reads
        count = 0;
        while (iter.hasNext()) {
            GATKRead rec = iter.next();
            assertTrue(rec.isUnmapped());
            count++;
        }
        assertEquals(count, 1000);
    }

    @Test
    public void testCreateArtificialUnmappedRead() {
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader();
        final GATKRead unmappedRead = ArtificialReadUtils.createArtificialUnmappedRead(header, new byte[]{'A'}, new byte[]{30});

        Assert.assertTrue(unmappedRead.isUnmapped());
        Assert.assertEquals(unmappedRead.getAssignedContig(), ReadConstants.UNSET_CONTIG);
        Assert.assertEquals(unmappedRead.getAssignedStart(), ReadConstants.UNSET_POSITION);
        Assert.assertEquals(unmappedRead.getBases(), new byte[]{'A'});
        Assert.assertEquals(unmappedRead.getBaseQualities(), new byte[]{30});
    }

    @Test
    public void testCreateArtificialUnmappedReadWithAssignedPosition() {
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader();
        final GATKRead unmappedRead = ArtificialReadUtils.createArtificialUnmappedReadWithAssignedPosition(header, "1", 50, new byte[]{'A'}, new byte[]{30});

        Assert.assertTrue(unmappedRead.isUnmapped());
        Assert.assertEquals(unmappedRead.getAssignedContig(), "1");
        Assert.assertEquals(unmappedRead.getAssignedStart(), 50);
        Assert.assertEquals(unmappedRead.getBases(), new byte[]{'A'});
        Assert.assertEquals(unmappedRead.getBaseQualities(), new byte[]{30});
    }

    @Test(dataProvider = "createNonIndelPileupElement")
    public void testCreateNonIndelPileupElement(final int offsetInRead, final Allele newAllele, final int lengthOfRead, final String gtBases) {
        final PileupElement pileupElement = ArtificialReadUtils.createNonIndelPileupElement(offsetInRead, newAllele, lengthOfRead);
        Assert.assertEquals(pileupElement.getRead().getBasesString(), gtBases);
    }

    @DataProvider(name="createNonIndelPileupElement")
    public Object[][] createNonIndelPileupElement() {
        return new Object[][] {
                {1, Allele.create("TTTT", true), 5, "ATTTT"},
                {0, Allele.create("TTTT", true), 5, "TTTTA"},
                {4, Allele.create("TTTT", true), 5, "AAAAT"},
                {1, Allele.create("TTTT", false), 5, "ATTTT"},
                {0, Allele.create("TTTT", false), 5, "TTTTA"},
                {4, Allele.create("TTTT", false), 5, "AAAAT"},
                {4, Allele.create("TCCT", false), 5, "AAAAT"},
                {4, Allele.create("CCCT", false), 5, "AAAAC"},
                {3, Allele.create("CCCT", false), 5, "AAACC"},
        };
    }

    @Test(dataProvider = "createInsertionPileupElement")
    public void testCreateInsertionPileup(final int offsetInRead, final Allele insertionAllele, final int lengthOfRead, final String gtBases) {
        final PileupElement pileupElement = ArtificialReadUtils.createSplicedInsertionPileupElement(offsetInRead, insertionAllele, lengthOfRead);
        Assert.assertEquals(pileupElement.getRead().getBasesString(), gtBases);
    }

    @DataProvider(name="createInsertionPileupElement")
    public Object[][] createInsertionileupElement() {
        return new Object[][] {
                {1, Allele.create("ATTT", true), 5, "AATTTAAA"},
                {0, Allele.create("ATTT", true), 5, "ATTTAAAA"},
                {4, Allele.create("ATTT", true), 5, "AAAAATTT"},
                {4, Allele.create("CTTT", true), 5, "AAAACTTT"},
                {1, Allele.create("ATTT", false), 5, "AATTTAAA"},
                {0, Allele.create("ATTT", false), 5, "ATTTAAAA"},
                {4, Allele.create("ATTT", false), 5, "AAAAATTT"},
        };
    }

    @Test(dataProvider = "createDeletionPileupElement")
    public void testCreateDeletionPileup(final int offsetInRead, final Allele referenceAllele, final int lengthOfRead, final String gtBases) {
        final PileupElement pileupElement = ArtificialReadUtils.createSplicedDeletionPileupElement(offsetInRead, referenceAllele, lengthOfRead);
        Assert.assertEquals(pileupElement.getRead().getBasesString(), gtBases);
    }

    @DataProvider(name="createDeletionPileupElement")
    public Object[][] createDeletionileupElement() {
        return new Object[][] {
                {1, Allele.create("AAAA", true), 5, "AA"},
                {0, Allele.create("CAAA", true), 5, "CA"},
                {1, Allele.create("AAAA", false), 5, "AA"},
                {0, Allele.create("CAAA", false), 5, "CA"},
                {4, Allele.create("AAAA", false), 5, "AAAAA"},
                {4, Allele.create("CAAA", false), 5, "AAAAC"}
        };
    }
}
