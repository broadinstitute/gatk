package org.broadinstitute.hellbender.utils;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

/**
 * Basic unit test for GenomeLoc
 */
public final class GenomeLocUnitTest extends GATKBaseTest {

    /**
     * Tests that we got a string parameter in correctly
     */
    @Test
    public void testIsBetween() {
        GenomeLoc locMiddle = hg19GenomeLocParser.createGenomeLoc("1", 3, 3);

        GenomeLoc locLeft = hg19GenomeLocParser.createGenomeLoc("1", 1, 1);
        GenomeLoc locRight = hg19GenomeLocParser.createGenomeLoc("1", 5, 5);

        Assert.assertTrue(locMiddle.isBetween(locLeft, locRight));
        Assert.assertFalse(locLeft.isBetween(locMiddle, locRight));
        Assert.assertFalse(locRight.isBetween(locLeft, locMiddle));

    }
    @Test
    public void testContigIndex() {
        GenomeLoc locOne = hg19GenomeLocParser.createGenomeLoc("1",1,1);
        Assert.assertEquals(0, locOne.getContigIndex());
        Assert.assertEquals("1", locOne.getContig());

        GenomeLoc locFour = hg19GenomeLocParser.createGenomeLoc("4",1,1);
        Assert.assertEquals(3, locFour.getContigIndex());
        Assert.assertEquals("4", locFour.contigName);

        GenomeLoc locNumber = hg19GenomeLocParser.createGenomeLoc(hg19Header.getSequenceDictionary().getSequence(0).getSequenceName(), 1, 1);
        Assert.assertEquals(0, locNumber.getContigIndex());
        Assert.assertEquals("1", locNumber.getContig());
        Assert.assertEquals(0, locOne.compareTo(locNumber));

    }

    @Test
    public void testCompareTo() {
        GenomeLoc twoOne = hg19GenomeLocParser.createGenomeLoc("2", 1);
        GenomeLoc twoFive = hg19GenomeLocParser.createGenomeLoc("2", 5);
        GenomeLoc twoOtherFive = hg19GenomeLocParser.createGenomeLoc("2", 5);
        Assert.assertEquals(twoFive.compareTo(twoOtherFive), 0);

        Assert.assertEquals(twoOne.compareTo(twoFive), -1);
        Assert.assertEquals(twoFive.compareTo(twoOne), 1);

        GenomeLoc oneOne = hg19GenomeLocParser.createGenomeLoc("1", 5);
        Assert.assertEquals(oneOne.compareTo(twoOne), -1);
        Assert.assertEquals(twoOne.compareTo(oneOne), 1);
    }

    @Test
    public void testEndpointSpan() throws Exception {
        GenomeLoc twoOne = hg19GenomeLocParser.createGenomeLoc("2", 1);
        GenomeLoc twoFive = hg19GenomeLocParser.createGenomeLoc("2", 5);
        final GenomeLoc twoOne_twoFive = twoOne.endpointSpan(twoFive);
        final GenomeLoc twoFive_twoOne = twoFive.endpointSpan(twoOne);

        final GenomeLoc expectedSpan =  hg19GenomeLocParser.createGenomeLoc("2", 1, 5);
        Assert.assertEquals(twoOne_twoFive, expectedSpan);
        Assert.assertEquals(twoFive_twoOne, expectedSpan);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testEndpointSpanUnmapped() throws Exception {
        GenomeLoc twoOne = hg19GenomeLocParser.createGenomeLoc("2", 1);
        GenomeLoc unmapped = GenomeLoc.UNMAPPED;
        twoOne.endpointSpan(unmapped);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testEndpointWrongContig() throws Exception {
        GenomeLoc twoOne = hg19GenomeLocParser.createGenomeLoc("2", 1);
        GenomeLoc threeFive = hg19GenomeLocParser.createGenomeLoc("3", 5);
        twoOne.endpointSpan(threeFive);
    }

    @Test
    public void testUnmappedSort() {
        GenomeLoc chr1 = hg19GenomeLocParser.createGenomeLoc("1",1,10000000);
        GenomeLoc chr2 = hg19GenomeLocParser.createGenomeLoc("2",1,10000000);
        GenomeLoc unmapped = GenomeLoc.UNMAPPED;

        List<GenomeLoc> unmappedOnly = Arrays.asList(unmapped);
        Collections.sort(unmappedOnly);
        Assert.assertEquals(unmappedOnly.size(),1,"Wrong number of elements in unmapped-only list.");
        Assert.assertEquals(unmappedOnly.get(0),unmapped,"List sorted in wrong order");

        List<GenomeLoc> chr1Presorted = Arrays.asList(chr1,unmapped);
        Collections.sort(chr1Presorted);
        Assert.assertEquals(chr1Presorted.size(),2,"Wrong number of elements in chr1,unmapped list.");
        Assert.assertEquals(chr1Presorted,Arrays.asList(chr1,unmapped),"List sorted in wrong order");

        List<GenomeLoc> chr1Inverted = Arrays.asList(unmapped,chr1);
        Collections.sort(chr1Inverted);
        Assert.assertEquals(chr1Inverted.size(),2,"Wrong number of elements in chr1,unmapped list.");
        Assert.assertEquals(chr1Inverted,Arrays.asList(chr1,unmapped),"List sorted in wrong order");

        List<GenomeLoc> chr1and2Presorted = Arrays.asList(chr1,chr2,unmapped);
        Collections.sort(chr1and2Presorted);
        Assert.assertEquals(chr1and2Presorted.size(),3,"Wrong number of elements in chr1,chr2,unmapped list.");
        Assert.assertEquals(chr1and2Presorted,Arrays.asList(chr1,chr2,unmapped),"List sorted in wrong order");

        List<GenomeLoc> chr1and2UnmappedInFront = Arrays.asList(unmapped,chr1,chr2);
        Collections.sort(chr1and2UnmappedInFront);
        Assert.assertEquals(chr1and2UnmappedInFront.size(),3,"Wrong number of elements in unmapped,chr1,chr2 list.");
        Assert.assertEquals(chr1and2UnmappedInFront,Arrays.asList(chr1,chr2,unmapped),"List sorted in wrong order");

        List<GenomeLoc> chr1and2UnmappedSandwiched = Arrays.asList(chr1,unmapped,chr2);
        Collections.sort(chr1and2UnmappedSandwiched);
        Assert.assertEquals(chr1and2UnmappedSandwiched.size(),3,"Wrong number of elements in chr1,unmapped,chr2 list.");
        Assert.assertEquals(chr1and2UnmappedSandwiched,Arrays.asList(chr1,chr2,unmapped),"List sorted in wrong order");
    }

    @Test
    public void testUnmappedMerge() {
        GenomeLoc chr1 = hg19GenomeLocParser.createGenomeLoc("1",1,10000000);
        GenomeLoc unmapped = GenomeLoc.UNMAPPED;

        List<GenomeLoc> oneUnmappedOnly = Arrays.asList(unmapped);
        oneUnmappedOnly = IntervalUtils.sortAndMergeIntervals(hg19GenomeLocParser,oneUnmappedOnly, IntervalMergingRule.OVERLAPPING_ONLY).toList();
        Assert.assertEquals(oneUnmappedOnly.size(),1,"Wrong number of elements in list.");
        Assert.assertEquals(oneUnmappedOnly.get(0),unmapped,"List sorted in wrong order");

        List<GenomeLoc> twoUnmapped = Arrays.asList(unmapped,unmapped);
        twoUnmapped = IntervalUtils.sortAndMergeIntervals(hg19GenomeLocParser,twoUnmapped,IntervalMergingRule.OVERLAPPING_ONLY).toList();
        Assert.assertEquals(twoUnmapped.size(),1,"Wrong number of elements in list.");
        Assert.assertEquals(twoUnmapped.get(0),unmapped,"List sorted in wrong order");

        List<GenomeLoc> twoUnmappedAtEnd = Arrays.asList(chr1,unmapped,unmapped);
        twoUnmappedAtEnd = IntervalUtils.sortAndMergeIntervals(hg19GenomeLocParser,twoUnmappedAtEnd,IntervalMergingRule.OVERLAPPING_ONLY).toList();
        Assert.assertEquals(twoUnmappedAtEnd.size(),2,"Wrong number of elements in list.");
        Assert.assertEquals(twoUnmappedAtEnd,Arrays.asList(chr1,unmapped),"List sorted in wrong order");

        List<GenomeLoc> twoUnmappedMixed = Arrays.asList(unmapped,chr1,unmapped);
        twoUnmappedMixed = IntervalUtils.sortAndMergeIntervals(hg19GenomeLocParser,twoUnmappedMixed,IntervalMergingRule.OVERLAPPING_ONLY).toList();
        Assert.assertEquals(twoUnmappedMixed.size(),2,"Wrong number of elements in list.");
        Assert.assertEquals(twoUnmappedMixed,Arrays.asList(chr1,unmapped),"List sorted in wrong order");
    }

    // -------------------------------------------------------------------------------------
    //
    // testing overlap detection
    //
    // -------------------------------------------------------------------------------------

    private class ReciprocalOverlapProvider extends TestDataProvider {
        GenomeLoc gl1, gl2;
        int overlapSize;
        double overlapFraction;

        private ReciprocalOverlapProvider(int start1, int stop1, int start2, int stop2) {
            super(ReciprocalOverlapProvider.class);
            gl1 = hg19GenomeLocParser.createGenomeLoc("1", start1, stop1);
            gl2 = hg19GenomeLocParser.createGenomeLoc("1", start2, stop2);

            int shared = 0;
            for ( int i = start1; i <= stop1; i++ ) {
                if ( i >= start2 && i <= stop2 )
                    shared++;
            }

            this.overlapSize = shared;
            this.overlapFraction = Math.min((1.0*shared)/gl1.size(), (1.0*shared)/gl2.size());
            super.setName(String.format("%d-%d / %d-%d overlap=%d / %.2f", start1, stop1, start2, stop2, overlapSize, overlapFraction));
        }
    }

    @DataProvider(name = "ReciprocalOverlapProvider")
    public Object[][] makeReciprocalOverlapProvider() {
        for ( int start1 = 1; start1 <= 10; start1++ ) {
            for ( int stop1 = start1; stop1 <= 10; stop1++ ) {
                new ReciprocalOverlapProvider(start1, stop1, 1, 10);
                new ReciprocalOverlapProvider(start1, stop1, 5, 10);
                new ReciprocalOverlapProvider(start1, stop1, 5, 7);
                new ReciprocalOverlapProvider(start1, stop1, 5, 15);
                new ReciprocalOverlapProvider(start1, stop1, 11, 20);

                new ReciprocalOverlapProvider(1, 10, start1, stop1);
                new ReciprocalOverlapProvider(5, 10, start1, stop1);
                new ReciprocalOverlapProvider(5, 7, start1, stop1);
                new ReciprocalOverlapProvider(5, 15, start1, stop1);
                new ReciprocalOverlapProvider(11, 20, start1, stop1);
            }
        }

        return ReciprocalOverlapProvider.getTests(ReciprocalOverlapProvider.class);
    }

    @Test(dataProvider = "ReciprocalOverlapProvider")
    public void testReciprocalOverlapProvider(ReciprocalOverlapProvider cfg) {
        if ( cfg.overlapSize == 0 ) {
            Assert.assertFalse(cfg.gl1.overlapsP(cfg.gl2));
        } else {
            Assert.assertTrue(cfg.gl1.overlapsP(cfg.gl2));
            Assert.assertEquals(cfg.gl1.intersect(cfg.gl2).size(), cfg.overlapSize);
            Assert.assertEquals(cfg.gl1.reciprocialOverlapFraction(cfg.gl2), cfg.overlapFraction);
        }
    }

    // -------------------------------------------------------------------------------------
    //
    // testing comparison, hashcode, and equals
    //
    // -------------------------------------------------------------------------------------

    @DataProvider(name = "GenomeLocComparisons")
    public Object[][] createGenomeLocComparisons() {
        List<Object[]> tests = new ArrayList<>();

        final int start = 10;
        for ( int stop = start; stop < start + 3; stop++ ) {
            final GenomeLoc g1 = hg19GenomeLocParser.createGenomeLoc("2", start, stop);
            for ( final String contig : Arrays.asList("1", "2", "3")) {
                for ( int start2 = start - 1; start2 <= stop + 1; start2++ ) {
                    for ( int stop2 = start2; stop2 < stop + 2; stop2++ ) {
                        final GenomeLoc g2 = hg19GenomeLocParser.createGenomeLoc(contig, start2, stop2);

                        ComparisonResult cmp = ComparisonResult.EQUALS;
                        if ( contig.equals("3") ) cmp = ComparisonResult.LESS_THAN;
                        else if ( contig.equals("1") ) cmp = ComparisonResult.GREATER_THAN;
                        else if ( start < start2 ) cmp = ComparisonResult.LESS_THAN;
                        else if ( start > start2 ) cmp = ComparisonResult.GREATER_THAN;
                        else if ( stop < stop2 ) cmp = ComparisonResult.LESS_THAN;
                        else if ( stop > stop2 ) cmp = ComparisonResult.GREATER_THAN;

                        tests.add(new Object[]{g1, g2, cmp});
                    }
                }
            }
        }

        return tests.toArray(new Object[][]{});
    }

    private enum ComparisonResult {
        LESS_THAN(-1),
        EQUALS(0),
        GREATER_THAN(1);

        final int cmp;

        private ComparisonResult(int cmp) {
            this.cmp = cmp;
        }
    }

    @Test(dataProvider = "GenomeLocComparisons")
    public void testGenomeLocComparisons(GenomeLoc g1, GenomeLoc g2, ComparisonResult expected) {
        Assert.assertEquals(g1.compareTo(g2), expected.cmp, "Comparing genome locs failed");
        Assert.assertEquals(g1.equals(g2), expected == ComparisonResult.EQUALS);
        if ( expected == ComparisonResult.EQUALS )
            Assert.assertEquals(g1.hashCode(), g2.hashCode(), "Equal genome locs don't have the same hash code");
    }

    // -------------------------------------------------------------------------------------
    //
    // testing merging functionality
    //
    // -------------------------------------------------------------------------------------

    private static final GenomeLoc loc1 = new GenomeLoc("1", 0, 10, 20);
    private static final GenomeLoc loc2 = new GenomeLoc("1", 0, 21, 30);
    private static final GenomeLoc loc3 = new GenomeLoc("1", 0, 31, 40);

    private class MergeTest {
        public List<GenomeLoc> locs;

        private MergeTest(final List<GenomeLoc> locs) {
            this.locs = locs;
        }
    }

    @DataProvider(name = "SGLtest")
    public Object[][] createFindVariantRegionsData() {
        List<Object[]> tests = new ArrayList<>();

        tests.add(new Object[]{new MergeTest(Arrays.<GenomeLoc>asList(loc1))});
        tests.add(new Object[]{new MergeTest(Arrays.<GenomeLoc>asList(loc1, loc2))});
        tests.add(new Object[]{new MergeTest(Arrays.<GenomeLoc>asList(loc1, loc2, loc3))});

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "SGLtest")
    public void testSimpleGenomeLoc(MergeTest test) {
        testMerge(test.locs);
    }

    @Test(expectedExceptions = GATKException.class)
    public void testNotContiguousLocs() {
        final List<GenomeLoc> locs = new ArrayList<>(1);
        locs.add(loc1);
        locs.add(loc3);
        testMerge(locs);
    }

    private void testMerge(final List<GenomeLoc> locs) {
        GenomeLoc result1 = locs.get(0);
        for ( int i = 1; i < locs.size(); i++ )
            result1 = GenomeLoc.merge(result1, locs.get(i));

        GenomeLoc result2 = GenomeLoc.merge(new TreeSet<>(locs));
        Assert.assertEquals(result1, result2);
        Assert.assertEquals(result1.getStart(), locs.get(0).getStart());
        Assert.assertEquals(result1.getStop(), locs.get(locs.size() - 1).getStop());
    }

    // -------------------------------------------------------------------------------------
    //
    // testing distance functionality
    //
    // -------------------------------------------------------------------------------------

    @Test
    public void testDistanceAcrossContigs() {
        final int chrSize = 1000;
        SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader(10, 0, chrSize);
        GenomeLocParser parser = new GenomeLocParser(header.getSequenceDictionary());
        GenomeLoc loc1 = parser.createGenomeLoc("3", 500);  // to check regular case
        GenomeLoc loc2 = parser.createGenomeLoc("7", 200);  // to check regular case
        GenomeLoc loc3 = parser.createGenomeLoc("0", 1);    // to check corner case
        GenomeLoc loc4 = parser.createGenomeLoc("9", 1000);// to check corner case
        GenomeLoc loc5 = parser.createGenomeLoc("7", 500);  // to make sure it does the right thing when in the same chromosome

        GenomeLoc loc6 = parser.createGenomeLoc("7", 200, 300);
        GenomeLoc loc7 = parser.createGenomeLoc("7", 500, 600);
        GenomeLoc loc8 = parser.createGenomeLoc("9", 500, 600);

        // Locus comparisons
        Assert.assertEquals(loc1.distanceAcrossContigs(loc2, header), 3*chrSize + chrSize-loc1.getStop() + loc2.getStart()); // simple case, smaller first
        Assert.assertEquals(loc2.distanceAcrossContigs(loc1, header), 3*chrSize + chrSize-loc1.getStop() + loc2.getStart()); // simple case, bigger first

        Assert.assertEquals(loc3.distanceAcrossContigs(loc4, header), 10*chrSize - 1); // corner case, smaller first
        Assert.assertEquals(loc4.distanceAcrossContigs(loc3, header), 10*chrSize - 1); // corner case, bigger first

        Assert.assertEquals(loc2.distanceAcrossContigs(loc5, header), 300); // same contig, smaller first
        Assert.assertEquals(loc5.distanceAcrossContigs(loc2, header), 300); // same contig, bigger first

        // Interval comparisons
        Assert.assertEquals(loc6.distanceAcrossContigs(loc7, header), 200); // same contig, smaller first
        Assert.assertEquals(loc7.distanceAcrossContigs(loc6, header), 200); // same contig, bigger first

        Assert.assertEquals(loc7.distanceAcrossContigs(loc8, header), chrSize + chrSize-loc7.stop + loc8.getStart()); // across contigs, smaller first
        Assert.assertEquals(loc8.distanceAcrossContigs(loc7, header), chrSize + chrSize-loc7.stop + loc8.getStart()); // across congits, bigger first
    }

}
