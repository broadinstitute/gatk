package org.broadinstitute.hellbender.utils;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.*;

import static org.testng.Assert.*;

public final class GenomeLocSortedSetUnitTest extends GATKBaseTest {

    private GenomeLocSortedSet mSortedSet = null;
    private SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader(NUMBER_OF_CHROMOSOMES, STARTING_CHROMOSOME, CHROMOSOME_SIZE);
    private static final int NUMBER_OF_CHROMOSOMES = 5;
    private static final int STARTING_CHROMOSOME = 1;
    private static final int CHROMOSOME_SIZE = 1000;

    private GenomeLocParser genomeLocParser;
    private String contigOneName;

    @BeforeClass
    public void setup() {
        genomeLocParser = new GenomeLocParser(header.getSequenceDictionary());
        contigOneName = header.getSequenceDictionary().getSequence(1).getSequenceName();
    }

    @BeforeMethod
    public void initializeSortedSet() {
        mSortedSet = new GenomeLocSortedSet(genomeLocParser);
    }

    @Test
    public void testAdd() {
        GenomeLoc g = genomeLocParser.createGenomeLoc(contigOneName, 0, 0);
        assertTrue(mSortedSet.size() == 0);
        mSortedSet.add(g);
        assertTrue(mSortedSet.size() == 1);
    }

    @Test
    public void testRemove() {
        assertTrue(mSortedSet.size() == 0);
        GenomeLoc g = genomeLocParser.createGenomeLoc(contigOneName, 0, 0);
        mSortedSet.add(g);
        assertTrue(mSortedSet.size() == 1);
        mSortedSet.remove(g);
        assertTrue(mSortedSet.size() == 0);
    }

    @Test
    public void addRegion() {
        assertTrue(mSortedSet.size() == 0);
        GenomeLoc g = genomeLocParser.createGenomeLoc(contigOneName, 1, 50);
        mSortedSet.add(g);
        GenomeLoc f = genomeLocParser.createGenomeLoc(contigOneName, 30, 80);
        mSortedSet.addRegion(f);
        assertTrue(mSortedSet.size() == 1);
    }

    @Test
    public void addRegionsOutOfOrder() {
        final String contigTwoName = header.getSequenceDictionary().getSequence(2).getSequenceName();
        assertTrue(mSortedSet.size() == 0);
        GenomeLoc g = genomeLocParser.createGenomeLoc(contigTwoName, 1, 50);
        mSortedSet.add(g);
        GenomeLoc f = genomeLocParser.createGenomeLoc(contigOneName, 30, 80);
        mSortedSet.addRegion(f);
        assertTrue(mSortedSet.size() == 2);
        assertTrue(mSortedSet.toList().get(0).getContig().equals(contigOneName));
        assertTrue(mSortedSet.toList().get(1).getContig().equals(contigTwoName));
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void addThrowsException() {
        assertTrue(mSortedSet.size() == 0);
        GenomeLoc g = genomeLocParser.createGenomeLoc(contigOneName, 1, 50);
        mSortedSet.add(g);
        GenomeLoc f = genomeLocParser.createGenomeLoc(contigOneName, 30, 80);
        mSortedSet.add(f);
    }

    @Test(expectedExceptions=IllegalArgumentException.class)
    public void testAddDuplicate() {
        assertTrue(mSortedSet.size() == 0);
        GenomeLoc g = genomeLocParser.createGenomeLoc(contigOneName, 0, 0);
        mSortedSet.add(g);
        assertTrue(mSortedSet.size() == 1);
        mSortedSet.add(g);
    }

    @Test
    public void mergingOverlappingBelow() {
        GenomeLoc g = genomeLocParser.createGenomeLoc(contigOneName, 0, 50);
        GenomeLoc e = genomeLocParser.createGenomeLoc(contigOneName, 49, 100);
        assertTrue(mSortedSet.size() == 0);
        mSortedSet.add(g);
        assertTrue(mSortedSet.size() == 1);
        mSortedSet.addRegion(e);
        assertTrue(mSortedSet.size() == 1);
        Iterator<GenomeLoc> iter = mSortedSet.iterator();
        GenomeLoc loc = iter.next();
        assertEquals(loc.getStart(), 0);
        assertEquals(loc.getStop(), 100);
        assertEquals(loc.getContigIndex(), 1);
    }

    @Test
    public void overlap() {
        for ( int i = 1; i < 6; i++ ) {
            final int start = i * 10;
            mSortedSet.add(genomeLocParser.createGenomeLoc(contigOneName, start, start + 1));
        }

        // test matches in and around interval
        assertFalse(mSortedSet.overlaps(genomeLocParser.createGenomeLoc(contigOneName, 9, 9)));
        assertTrue(mSortedSet.overlaps(genomeLocParser.createGenomeLoc(contigOneName, 10, 10)));
        assertTrue(mSortedSet.overlaps(genomeLocParser.createGenomeLoc(contigOneName, 11, 11)));
        assertFalse(mSortedSet.overlaps(genomeLocParser.createGenomeLoc(contigOneName, 12, 12)));

        // test matches spanning intervals
        assertTrue(mSortedSet.overlaps(genomeLocParser.createGenomeLoc(contigOneName, 14, 20)));
        assertTrue(mSortedSet.overlaps(genomeLocParser.createGenomeLoc(contigOneName, 11, 15)));
        assertTrue(mSortedSet.overlaps(genomeLocParser.createGenomeLoc(contigOneName, 30, 40)));
        assertTrue(mSortedSet.overlaps(genomeLocParser.createGenomeLoc(contigOneName, 51, 53)));

        // test miss
        assertFalse(mSortedSet.overlaps(genomeLocParser.createGenomeLoc(contigOneName, 12, 19)));

        // test exact match after miss
        assertTrue(mSortedSet.overlaps(genomeLocParser.createGenomeLoc(contigOneName, 40, 41)));

        // test matches at beginning of intervals
        assertFalse(mSortedSet.overlaps(genomeLocParser.createGenomeLoc(contigOneName, 5, 6)));
        assertTrue(mSortedSet.overlaps(genomeLocParser.createGenomeLoc(contigOneName, 0, 10)));

        // test matches at end of intervals
        assertFalse(mSortedSet.overlaps(genomeLocParser.createGenomeLoc(contigOneName, 52, 53)));
        assertTrue(mSortedSet.overlaps(genomeLocParser.createGenomeLoc(contigOneName, 51, 53)));
        assertFalse(mSortedSet.overlaps(genomeLocParser.createGenomeLoc(contigOneName, 52, 53)));
    }

    @Test
    public void mergingOverlappingAbove() {
        GenomeLoc e = genomeLocParser.createGenomeLoc(contigOneName, 0, 50);
        GenomeLoc g = genomeLocParser.createGenomeLoc(contigOneName, 49, 100);
        assertTrue(mSortedSet.size() == 0);
        mSortedSet.add(g);
        assertTrue(mSortedSet.size() == 1);
        mSortedSet.addRegion(e);
        assertTrue(mSortedSet.size() == 1);
        Iterator<GenomeLoc> iter = mSortedSet.iterator();
        GenomeLoc loc = iter.next();
        assertEquals(loc.getStart(), 0);
        assertEquals(loc.getStop(), 100);
        assertEquals(loc.getContigIndex(), 1);
    }

    @Test
    public void deleteAllByRegion() {
        GenomeLoc e = genomeLocParser.createGenomeLoc(contigOneName, 1, 100);
        mSortedSet.add(e);
        for (int x = 1; x < 101; x++) {
            GenomeLoc del = genomeLocParser.createGenomeLoc(contigOneName,x,x);
            mSortedSet = mSortedSet.subtractRegions(new GenomeLocSortedSet(genomeLocParser,del));
        }
        assertTrue(mSortedSet.isEmpty());
    }

    @Test
    public void deleteSomeByRegion() {
        GenomeLoc e = genomeLocParser.createGenomeLoc(contigOneName, 1, 100);
        mSortedSet.add(e);
        for (int x = 1; x < 50; x++) {
            GenomeLoc del = genomeLocParser.createGenomeLoc(contigOneName,x,x);
            mSortedSet = mSortedSet.subtractRegions(new GenomeLocSortedSet(genomeLocParser,del));
        }
        assertTrue(!mSortedSet.isEmpty());
        assertTrue(mSortedSet.size() == 1);
        GenomeLoc loc = mSortedSet.iterator().next();
        assertTrue(loc.getStop() == 100);
        assertTrue(loc.getStart() == 50);

    }

    @Test
    public void deleteSuperRegion() {
        GenomeLoc e = genomeLocParser.createGenomeLoc(contigOneName, 10, 20);
        GenomeLoc g = genomeLocParser.createGenomeLoc(contigOneName, 70, 100);
        mSortedSet.add(g);
        mSortedSet.addRegion(e);
        assertTrue(mSortedSet.size() == 2);
        // now delete a region
        GenomeLoc d = genomeLocParser.createGenomeLoc(contigOneName, 15, 75);
        mSortedSet = mSortedSet.subtractRegions(new GenomeLocSortedSet(genomeLocParser,d));
        Iterator<GenomeLoc> iter = mSortedSet.iterator();
        GenomeLoc loc = iter.next();
        assertTrue(loc.getStart() == 10);
        assertTrue(loc.getStop() == 14);
        assertTrue(loc.getContigIndex() == 1);

        loc = iter.next();
        assertTrue(loc.getStart() == 76);
        assertTrue(loc.getStop() == 100);
        assertTrue(loc.getContigIndex() == 1);
    }

    @Test
    public void substractComplexExample() {
        GenomeLoc e = genomeLocParser.createGenomeLoc(contigOneName, 1, 20);
        mSortedSet.add(e);

        GenomeLoc r1 = genomeLocParser.createGenomeLoc(contigOneName, 3, 5);
        GenomeLoc r2 = genomeLocParser.createGenomeLoc(contigOneName, 10, 12);
        GenomeLoc r3 = genomeLocParser.createGenomeLoc(contigOneName, 16, 18);
        GenomeLocSortedSet toExclude = new GenomeLocSortedSet(genomeLocParser,Arrays.asList(r1, r2, r3));

        GenomeLocSortedSet remaining = mSortedSet.subtractRegions(toExclude);
//        logger.debug("Initial   " + mSortedSet);
//        logger.debug("Exclude   " + toExclude);
//        logger.debug("Remaining " + remaining);

        assertEquals(mSortedSet.coveredSize(), 20);
        assertEquals(toExclude.coveredSize(), 9);
        assertEquals(remaining.coveredSize(), 11);

        Iterator<GenomeLoc> it = remaining.iterator();
        GenomeLoc p1 = it.next();
        GenomeLoc p2 = it.next();
        GenomeLoc p3 = it.next();
        GenomeLoc p4 = it.next();

        assertEquals(genomeLocParser.createGenomeLoc(contigOneName, 1, 2), p1);
        assertEquals(genomeLocParser.createGenomeLoc(contigOneName, 6, 9), p2);
        assertEquals(genomeLocParser.createGenomeLoc(contigOneName, 13, 15), p3);
        assertEquals(genomeLocParser.createGenomeLoc(contigOneName, 19, 20), p4);
    }

    private void testSizeBeforeLocX(int pos, int size) {
        GenomeLoc test = genomeLocParser.createGenomeLoc(contigOneName, pos, pos);
        assertEquals(mSortedSet.sizeBeforeLoc(test), size, String.format("X pos=%d size=%d", pos, size));
    }

    @Test
    public void testSizeBeforeLoc() {
        GenomeLoc r1 = genomeLocParser.createGenomeLoc(contigOneName, 3, 5);
        GenomeLoc r2 = genomeLocParser.createGenomeLoc(contigOneName, 10, 12);
        GenomeLoc r3 = genomeLocParser.createGenomeLoc(contigOneName, 16, 18);
        mSortedSet.addAll(Arrays.asList(r1,r2,r3));

        testSizeBeforeLocX(2, 0);
        testSizeBeforeLocX(3, 0);
        testSizeBeforeLocX(4, 1);
        testSizeBeforeLocX(5, 2);
        testSizeBeforeLocX(6, 3);

        testSizeBeforeLocX(10, 3);
        testSizeBeforeLocX(11, 4);
        testSizeBeforeLocX(12, 5);
        testSizeBeforeLocX(13, 6);
        testSizeBeforeLocX(15, 6);

        testSizeBeforeLocX(16, 6);
        testSizeBeforeLocX(17, 7);
        testSizeBeforeLocX(18, 8);
        testSizeBeforeLocX(19, 9);
        testSizeBeforeLocX(50, 9);
        testSizeBeforeLocX(50, (int)mSortedSet.coveredSize());
    }


    @Test
    public void fromSequenceDictionary() {
        mSortedSet = GenomeLocSortedSet.createSetFromSequenceDictionary(this.header.getSequenceDictionary());
        // we should have sequence
        assertTrue(mSortedSet.size() == GenomeLocSortedSetUnitTest.NUMBER_OF_CHROMOSOMES);
        int seqNumber = 0;
        for (GenomeLoc loc : mSortedSet) {
            assertTrue(loc.getStart() == 1);
            assertTrue(loc.getStop() == GenomeLocSortedSetUnitTest.CHROMOSOME_SIZE);
            assertTrue(loc.getContigIndex() == seqNumber);
            ++seqNumber;
        }
        assertTrue(seqNumber == GenomeLocSortedSetUnitTest.NUMBER_OF_CHROMOSOMES);
    }

    // -----------------------------------------------------------------------------------------------
    //
    // Test getOverlapping
    //
    // -----------------------------------------------------------------------------------------------

    @DataProvider(name = "GetOverlapping")
    public Object[][] makeGetOverlappingTest() throws Exception {
        final GenomeLocParser genomeLocParser = new GenomeLocParser(new CachingIndexedFastaSequenceFile(
                IOUtils.getPath(exampleReference)));

        List<Object[]> tests = new ArrayList<>();

        final GenomeLoc prev1 = genomeLocParser.createGenomeLoc("1", 1, 10);
        final GenomeLoc prev2 = genomeLocParser.createGenomeLoc("1", 20, 50);
        final GenomeLoc post1 = genomeLocParser.createGenomeLoc("3", 1, 10);
        final GenomeLoc post2 = genomeLocParser.createGenomeLoc("3", 20, 50);

        final int chr20Length = genomeLocParser.getSequenceDictionary().getSequence("2").getSequenceLength();
        for ( final int regionStart : Arrays.asList(1, 10, chr20Length - 10, chr20Length) ) {
            for ( final int regionSize : Arrays.asList(1, 10, 100) ) {
                final GenomeLoc region = genomeLocParser.createGenomeLocOnContig("2", regionStart, regionStart + regionSize);
                final GenomeLoc spanning = genomeLocParser.createGenomeLocOnContig("2", regionStart - 10, region.getStop() + 10);
                final GenomeLoc before_into = genomeLocParser.createGenomeLocOnContig("2", regionStart - 10, regionStart + 1);
                final GenomeLoc middle = genomeLocParser.createGenomeLocOnContig("2", regionStart + 1, regionStart + 2);
                final GenomeLoc middle_past = genomeLocParser.createGenomeLocOnContig("2", region.getStop()-1, region.getStop()+10);

                final List<GenomeLoc> potentials = new LinkedList<>();
                potentials.add(region);
                if ( spanning != null ) potentials.add(spanning);
                if ( before_into != null ) potentials.add(before_into);
                if ( middle != null ) potentials.add(middle);
                if ( middle_past != null ) potentials.add(middle_past);

                for ( final int n : Arrays.asList(1, 2, 3) ) {
                    for ( final List<GenomeLoc> regions : Utils.makePermutations(potentials, n, false) ) {
                        tests.add(new Object[]{new GenomeLocSortedSet(genomeLocParser, regions), region});
                        tests.add(new Object[]{new GenomeLocSortedSet(genomeLocParser, Utils.append(regions, prev1)), region});
                        tests.add(new Object[]{new GenomeLocSortedSet(genomeLocParser, Utils.append(regions, prev1, prev2)), region});
                        tests.add(new Object[]{new GenomeLocSortedSet(genomeLocParser, Utils.append(regions, post1)), region});
                        tests.add(new Object[]{new GenomeLocSortedSet(genomeLocParser, Utils.append(regions, post1, post2)), region});
                        tests.add(new Object[]{new GenomeLocSortedSet(genomeLocParser, Utils.append(regions, prev1, post1)), region});
                        tests.add(new Object[]{new GenomeLocSortedSet(genomeLocParser, Utils.append(regions, prev1, prev2, post1, post2)), region});
                    }
                }
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "GetOverlapping")
    public void testGetOverlapping(final GenomeLocSortedSet intervals, final GenomeLoc region) {
        final List<GenomeLoc> expectedOverlapping = intervals.getOverlappingFullSearch(region);
        final List<GenomeLoc> actualOverlapping = intervals.getOverlapping(region);
        Assert.assertEquals(actualOverlapping, expectedOverlapping);
        Assert.assertEquals(intervals.overlaps(region), ! expectedOverlapping.isEmpty(), "GenomeLocSortedSet.overlaps didn't return expected result");
    }
}
