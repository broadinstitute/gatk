package org.broadinstitute.hellbender.tools.exome;

import htsjdk.samtools.SAMSequenceRecord;
import org.broadinstitute.hellbender.utils.GenomeLoc;
import org.broadinstitute.hellbender.utils.IndexRange;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;

/**
 * Unit tests for {@link IntervalBackedExonCollection}.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class IntervalBackedExonCollectionUnitTest extends BaseTest {

    /**
     * Average exon size for randomly generated test data.
     */
    private final int averageExonSize = 100;

    /**
     * Std. dev. for exon size for randomly generated test data.
     */
    private final float sdExonSize = 20;

    /**
     * Average inter-exon gap in bp for randomly generated test data.
     */
    private final int averageExonIntergapSize = 100;

    /**
     * Std. dev. for inter-exon gap in bp for randomly generated test data.
     */
    private final float sdExonIntergapSize = 20;

    /**
     * Minimum inter-exon gap size for randomly generated test data.
     */
    private final int minimumExonIntergapSize = 20;

    /**
     * Minimum exon size for randomly generated test data.
     */
    private final int minimumExonSize = 1;

    /**
     * Reference to the test intervals.
     * <p>
     *     Initialized before all tests in {@link #setUp}.
     * </p>
     */
    private List<GenomeLoc> nonOverlappingExomeIntervals;

    @BeforeClass
    public void setUp() {
        nonOverlappingExomeIntervals = new ArrayList<>(100);
        final Random rdn = new Random(13);// some "random" but fix seed to make sure errors are deterministic.
        for (int i = 0; i < ExomeToolsTestUtils.REFERENCE_DICTIONARY.size(); i++) {
            int current = 0;
            final SAMSequenceRecord sequence = ExomeToolsTestUtils.REFERENCE_DICTIONARY.getSequence(i);
            while (current < sequence.getSequenceLength()) {
                int start = current + Math.max(minimumExonIntergapSize,
                        (int)Math.round(rdn.nextGaussian() * sdExonIntergapSize + averageExonIntergapSize));
                if (start >= sequence.getSequenceLength()) {
                    break;
                }
                int size = Math.max(minimumExonSize,
                        (int)Math.round(rdn.nextGaussian() * sdExonSize + averageExonSize));
                int stop = start + size - 1;
                if (stop >= sequence.getSequenceLength()) {
                    break;
                }
                nonOverlappingExomeIntervals.add(
                        ExomeToolsTestUtils.GENOME_LOC_FACTORY.createGenomeLoc(
                                sequence.getSequenceName(), start, stop));
                current = stop + 1;
            }
        }
    }

    @Test
    public void testCorrectInitialization() {
        new IntervalBackedExonCollection(nonOverlappingExomeIntervals);
    }

    @Test(dataProvider="correctRangeObjectInitializationData")
    public void testCorrectRangeObjectInitialization(final int from, final int to) {
        new IndexRange(from,to);
    }

    @Test
    public void testCorrectEmptyInitialization() {
        new IntervalBackedExonCollection(Collections.emptyList());
    }

    @Test(expectedExceptions = {IllegalArgumentException.class})
    public void testIncorrectInitialization1() {
        new IntervalBackedExonCollection(null);
    }

    @Test(expectedExceptions = {IllegalArgumentException.class})
    public void testIncorrectInitialization2() {
        final List<GenomeLoc> intervals = new ArrayList<>(nonOverlappingExomeIntervals);
        intervals.add(nonOverlappingExomeIntervals.size() >> 1,null);
        new IntervalBackedExonCollection(intervals);
    }

    @Test(expectedExceptions = {IllegalArgumentException.class})
    public void testIncorrectInitialization3() {
        final List<GenomeLoc> intervals = new ArrayList<>(nonOverlappingExomeIntervals);
        intervals.add(GenomeLoc.UNMAPPED);
        new IntervalBackedExonCollection(intervals);
    }

    @Test(expectedExceptions = {IllegalArgumentException.class})
    public void testIncorrectInitialization4() {
        final List<GenomeLoc> intervals = new ArrayList<>(nonOverlappingExomeIntervals);
        intervals.add(nonOverlappingExomeIntervals.size() >> 1,GenomeLoc.WHOLE_GENOME);
        new IntervalBackedExonCollection(intervals);
    }

    @Test(expectedExceptions = {IllegalArgumentException.class})
    public void testIncorrectInitialization5() {
        final List<GenomeLoc> intervals = new ArrayList<>(nonOverlappingExomeIntervals);
        intervals.add(1,intervals.get(0));
        new IntervalBackedExonCollection(intervals);
    }

    @Test(expectedExceptions = {IllegalArgumentException.class})
    public void testIncorrectInitialization6() {
        final List<GenomeLoc> intervals = new ArrayList<>(nonOverlappingExomeIntervals);
        intervals.add(1, intervals.get(0).getStopLocation());
        new IntervalBackedExonCollection(intervals);
    }

    @Test(expectedExceptions = {IllegalArgumentException.class})
    public void testIncorrectInitialization7() {
        final List<GenomeLoc> intervals = new ArrayList<>(nonOverlappingExomeIntervals);
        final GenomeLoc last = intervals.remove(intervals.size() - 1);
        intervals.add(2, last);
        @SuppressWarnings("unused")
        final IntervalBackedExonCollection exonDB = new IntervalBackedExonCollection(intervals);
    }


    @Test(dependsOnMethods={"testCorrectInitialization"})
    public void testExonCount() {
        final IntervalBackedExonCollection exonDB = new IntervalBackedExonCollection(nonOverlappingExomeIntervals);
        Assert.assertEquals(exonDB.exonCount(),nonOverlappingExomeIntervals.size());
    }

    @Test(dependsOnMethods={"testCorrectInitialization"})
    public void testExomeSize() {
        int size = 0;
        for (final GenomeLoc loc : nonOverlappingExomeIntervals) {
            size += loc.size();
        }
        final IntervalBackedExonCollection exonDB = new IntervalBackedExonCollection(nonOverlappingExomeIntervals);
        Assert.assertEquals(exonDB.exomeSize(),size);
    }

    @Test(dependsOnMethods = {"testCorrectInitialization"})
    public void testExons() {
        final IntervalBackedExonCollection exonDB = new IntervalBackedExonCollection(nonOverlappingExomeIntervals);
        final List<GenomeLoc> exons = exonDB.exons();
        Assert.assertNotNull(exons);
        Assert.assertEquals(exons.size(),nonOverlappingExomeIntervals.size());
        for (int i = 0; i < exons.size(); i++) {
            Assert.assertEquals(exons.get(i), nonOverlappingExomeIntervals.get(i));
        }
    }

    @Test(dependsOnMethods = {"testCorrectInitialization"})
    public void testExonsFullRange() {
        final IntervalBackedExonCollection exonDB = new IntervalBackedExonCollection(nonOverlappingExomeIntervals);
        final List<GenomeLoc> exons = exonDB.exons(0,exonDB.exonCount());
        Assert.assertNotNull(exons);
        Assert.assertEquals(exons.size(),nonOverlappingExomeIntervals.size());
        for (int i = 0; i < exons.size(); i++) {
            Assert.assertEquals(exons.get(i), nonOverlappingExomeIntervals.get(i));
        }
    }

    @Test(dependsOnMethods = {"testCorrectInitialization","testCorrectRangeObjectInitialization"})
    public void testExonFullRangeObject() {
        final IntervalBackedExonCollection exonDB = new IntervalBackedExonCollection(nonOverlappingExomeIntervals);
        final IndexRange range = new IndexRange(0,exonDB.exonCount());
        final List<GenomeLoc> exons = exonDB.exons(range);
        Assert.assertNotNull(exons);
        Assert.assertEquals(exons.size(),nonOverlappingExomeIntervals.size());
        for (int i = 0; i < exons.size(); i++) {
            Assert.assertEquals(exons.get(i), nonOverlappingExomeIntervals.get(i));
        }
    }

    @Test(dependsOnMethods = {"testCorrectInitialization"})
    public void testEmptyRange() {
        final IntervalBackedExonCollection exonDB = new IntervalBackedExonCollection(nonOverlappingExomeIntervals);
        final List<GenomeLoc> exons = exonDB.exons(0,0);
        Assert.assertNotNull(exons);
        Assert.assertEquals(exons.size(),0);
    }

    @Test(dependsOnMethods = {"testCorrectInitialization","testCorrectRangeObjectInitialization"})
    public void testEmptyRangeObject() {
        final IntervalBackedExonCollection exonDB = new IntervalBackedExonCollection(nonOverlappingExomeIntervals);
        final IndexRange range = new IndexRange(0,0);
        final List<GenomeLoc> exons = exonDB.exons(range);
        Assert.assertNotNull(exons);
        Assert.assertEquals(exons.size(),0);
    }

    @Test(dependsOnMethods = {"testCorrectInitialization"},
          dataProvider = "arbitraryExonRageData")
    public void testArbitraryExonRange(final int from, final int to) {
        final IntervalBackedExonCollection exonDB = new IntervalBackedExonCollection(nonOverlappingExomeIntervals);
        final List<GenomeLoc> exons = exonDB.exons(from,to);
        Assert.assertNotNull(exons);
        final int expectedSize = to - from;
        Assert.assertEquals(exons.size(),expectedSize);
        for (int i = from; i < to; i++) {
            Assert.assertEquals(exons.get(i - from), nonOverlappingExomeIntervals.get(i));
        }
    }

    @Test(dependsOnMethods = {"testCorrectInitialization","testCorrectRangeObjectInitialization"},
            dataProvider = "arbitraryExonRageData")
    public void testArbitraryExonRangeObject(final int from, final int to) {
        final IntervalBackedExonCollection exonDB = new IntervalBackedExonCollection(nonOverlappingExomeIntervals);
        final IndexRange range = new IndexRange(from,to);
        final List<GenomeLoc> exons = exonDB.exons(range);
        Assert.assertNotNull(exons);
        final int expectedSize = to - from;
        Assert.assertEquals(exons.size(), expectedSize);
        for (int i = from; i < to; i++) {
            Assert.assertEquals(exons.get(i - from), nonOverlappingExomeIntervals.get(i));
        }
    }

    @Test(dependsOnMethods = "testCorrectInitialization")
    public void testLocation() {
        final IntervalBackedExonCollection exonDB = new IntervalBackedExonCollection(nonOverlappingExomeIntervals);
        for (final GenomeLoc loc : nonOverlappingExomeIntervals) {
            Assert.assertEquals(exonDB.location(loc), loc);
        }
        Assert.assertEquals(exonDB.location(GenomeLoc.WHOLE_GENOME), GenomeLoc.WHOLE_GENOME);
        Assert.assertEquals(exonDB.location(GenomeLoc.UNMAPPED), GenomeLoc.UNMAPPED);
    }

    @Test(dependsOnMethods = "testCorrectInitialization",
        expectedExceptions = IllegalArgumentException.class )
    public void testWrongLocation() {
        final IntervalBackedExonCollection exonDB = new IntervalBackedExonCollection(nonOverlappingExomeIntervals);
        exonDB.location(null);
    }

    @Test(dependsOnMethods = "testCorrectInitialization")
    public void testLocationFromIndex() {
        final IntervalBackedExonCollection exonDB = new IntervalBackedExonCollection(nonOverlappingExomeIntervals);
        for (int i = 0; i < nonOverlappingExomeIntervals.size(); i++) {
            Assert.assertEquals(exonDB.location(i), nonOverlappingExomeIntervals.get(i));
        }
        Assert.assertEquals(exonDB.location(GenomeLoc.WHOLE_GENOME),GenomeLoc.WHOLE_GENOME);
        Assert.assertEquals(exonDB.location(GenomeLoc.UNMAPPED), GenomeLoc.UNMAPPED);
    }

    @Test(dependsOnMethods = {"testCorrectInitialization"},
          expectedExceptions = {IndexOutOfBoundsException.class, IllegalArgumentException.class},
          dataProvider = "invalidRangeData")
    public void testInvalidRange(final int from, final int to) {
        final IntervalBackedExonCollection exonDB = new IntervalBackedExonCollection(nonOverlappingExomeIntervals);
        exonDB.exons(from, to);
    }

    @Test(dependsOnMethods = {"testCorrectInitialization","testCorrectRangeObjectInitialization"},
            expectedExceptions = {IndexOutOfBoundsException.class, IllegalArgumentException.class},
            dataProvider = "invalidRangeData")
    public void testInvalidRangeUsingObject(final int from, final int to) {
        final IntervalBackedExonCollection exonDB = new IntervalBackedExonCollection(nonOverlappingExomeIntervals);
        final IndexRange range = new IndexRange(from,to);
        exonDB.exons(range);
    }

    @Test(dependsOnMethods = "testCorrectInitialization",
          dataProvider = "exonLookUpData")
    public void testExonByLocation(final GenomeLoc location, final GenomeLoc expected,
                                   @SuppressWarnings("unused") final int index) {
        final IntervalBackedExonCollection exonDB = new IntervalBackedExonCollection(nonOverlappingExomeIntervals);
        final GenomeLoc actual = exonDB.exon(location);
        if (expected == null) {
            Assert.assertNull(actual);
        } else {
            Assert.assertNotNull(actual);
            Assert.assertEquals(actual,expected);
        }
    }

    @Test(dependsOnMethods = "testCorrectInitialization",
              dataProvider = "exonLookUpData")
    public void testExonIndexByLocation(final GenomeLoc location, @SuppressWarnings("unused")  final GenomeLoc expected,
                                   final int index) {
        final IntervalBackedExonCollection exonDB = new IntervalBackedExonCollection(nonOverlappingExomeIntervals);
        final int actual = exonDB.index(location);
        Assert.assertEquals(actual, index);
    }

    @Test(dependsOnMethods = "testCorrectInitialization",
          expectedExceptions = ExonCollection.AmbiguousExonException.class)
    public void testIndexWholeGenomeWithManyExons() {
        final IntervalBackedExonCollection exonDB = new IntervalBackedExonCollection(nonOverlappingExomeIntervals);
        exonDB.index(GenomeLoc.WHOLE_GENOME);
    }

    @Test(dependsOnMethods = "testCorrectInitialization")
    public void testIndexWholeGenomeWithNoExon() {
        final IntervalBackedExonCollection exonDB = new IntervalBackedExonCollection(Collections.emptyList());
        Assert.assertEquals(exonDB.index(GenomeLoc.WHOLE_GENOME), -1);
    }

    @Test(dependsOnMethods = "testCorrectInitialization")
    public void testIndexWholeGenomeWithOneExon() {
        final IntervalBackedExonCollection exonDB = new IntervalBackedExonCollection(nonOverlappingExomeIntervals.subList(0, 1));
        Assert.assertEquals(exonDB.index(GenomeLoc.WHOLE_GENOME), 0);
    }

    @Test(dependsOnMethods = "testCorrectInitialization",
              dataProvider = "wrongExonLookUpData",
            expectedExceptions = {ExonCollection.AmbiguousExonException.class, IllegalArgumentException.class})
    public void testWrongExonIndexByLocation(final GenomeLoc location) {
        final IntervalBackedExonCollection exonDB = new IntervalBackedExonCollection(nonOverlappingExomeIntervals);
        exonDB.index(location);
    }

    @Test(dependsOnMethods = "testCorrectInitialization",
            dataProvider = "wrongExonLookUpData",
            expectedExceptions = {ExonCollection.AmbiguousExonException.class, IllegalArgumentException.class})
    public void testWrongExonByLocation(final GenomeLoc location) {
        final IntervalBackedExonCollection exonDB = new IntervalBackedExonCollection(nonOverlappingExomeIntervals);
        exonDB.exon(location);
    }

    @Test(dependsOnMethods = "testCorrectInitialization")
    public void testExonByName() {
        final IntervalBackedExonCollection exonDB = new IntervalBackedExonCollection(nonOverlappingExomeIntervals);

        for (final GenomeLoc loc : nonOverlappingExomeIntervals) {
            Assert.assertEquals(exonDB.exon(loc.toString()), loc);
        }

        Assert.assertNull(exonDB.exon("no-id"));
    }

    @Test(dependsOnMethods = "testCorrectInitialization")
    public void testExonIndexByName() {
        final IntervalBackedExonCollection exonDB = new IntervalBackedExonCollection(nonOverlappingExomeIntervals);

        for (int i = 0; i < nonOverlappingExomeIntervals.size(); i++) {
            Assert.assertEquals(exonDB.index(nonOverlappingExomeIntervals.get(i)),i);
        }

        Assert.assertEquals(exonDB.index("no-id"), -1);
    }

    @Test(dependsOnMethods = {"testCorrectInitialization"},
        dataProvider="indexRangeData")
    public void testIndexRange(final GenomeLoc region, final IndexRange expected) {
        final IntervalBackedExonCollection exonDB = new IntervalBackedExonCollection(nonOverlappingExomeIntervals);
        final IndexRange range = exonDB.indexRange(region);
        Assert.assertEquals(range, expected);
    }

    @Test(dependsOnMethods = {"testCorrectInitialization"},
            dataProvider="indexRangeData")
    public void testForEachExon(final GenomeLoc region, final IndexRange expected) {
        final IntervalBackedExonCollection exonDB = new IntervalBackedExonCollection(nonOverlappingExomeIntervals);
        final List<GenomeLoc> exons = exonDB.exons(region);
        final List<Integer> exonIndices = new ArrayList<>(exons.size());
        final List<GenomeLoc> exonObjects = new ArrayList<>(exons.size());
        exonDB.forEachExon(region,(index,exon) -> {
            exonIndices.add(index);
            exonObjects.add(exon);
        });
        Assert.assertEquals(exonIndices.size(),exons.size());
        for (int i = 0; i < exons.size(); i++) {
            Assert.assertEquals(exonIndices.get(i).intValue(),expected.from + i);
            Assert.assertEquals(exonObjects.get(i),exons.get(i));
        }
    }

    @Test(dependsOnMethods = {"testCorrectInitialization"},
       dataProvider = "indexRangeData")
    public void testExons(final GenomeLoc region, final IndexRange expected) {
        final IntervalBackedExonCollection exonDB = new IntervalBackedExonCollection(nonOverlappingExomeIntervals);
        final List<GenomeLoc> exons = exonDB.exons(region);
        final List<GenomeLoc> expectedExons =  nonOverlappingExomeIntervals.subList(expected.from,expected.to);
        Assert.assertEquals(exons,expectedExons);
    }

    @Test(dependsOnMethods={"testCorrectInitialization","testExonCount"})
    public void testExonName() {
        final IntervalBackedExonCollection exonDB = new IntervalBackedExonCollection(nonOverlappingExomeIntervals);
        for (int i = 0; i < exonDB.exonCount(); i++) {
            final GenomeLoc loc = exonDB.exon(i);
            final String name = exonDB.intervalName(loc);
            final GenomeLoc observed = exonDB.exon(name);
            Assert.assertEquals(observed,loc);
            Assert.assertEquals(exonDB.index(name),i);
        }
        Assert.assertNull(exonDB.exon("no-interval"));
        Assert.assertEquals(exonDB.index("no-interval"),-1);
    }

    @DataProvider(name="indexRangeData")
    public Object[][] indexRangeData() {
        final List<Object[]> result = new ArrayList<>();

        final Random rdn = new Random(131313);
        for (int i = 0; i < nonOverlappingExomeIntervals.size(); i++) {
            result.add(new Object[]{nonOverlappingExomeIntervals.get(i), new IndexRange(i, i + 1)});
            result.add(new Object[]{
                    ExomeToolsTestUtils.createGenomeLoc(nonOverlappingExomeIntervals.get(i).getContig(),
                            nonOverlappingExomeIntervals.get(i).getStart() - 1), new IndexRange(i,i)});
            int j;
            for (j = i + 1; j < nonOverlappingExomeIntervals.size() && rdn.nextBoolean(); j++) {
                if (!nonOverlappingExomeIntervals.get(j).getContig().equals(nonOverlappingExomeIntervals.get(i).getContig())) {
                    break;
                }
            }
            result.add(new Object[] {
                    ExomeToolsTestUtils.createGenomeLoc(nonOverlappingExomeIntervals.get(i).getContig(),
                            nonOverlappingExomeIntervals.get(i).getStart(),
                            nonOverlappingExomeIntervals.get(j - 1).getStop()),
                    new IndexRange(i, j)
            });
        }
        result.add(new Object[]{GenomeLoc.WHOLE_GENOME, new IndexRange(0, nonOverlappingExomeIntervals.size())});
        result.add(new Object[]{GenomeLoc.UNMAPPED, new IndexRange(nonOverlappingExomeIntervals.size(), nonOverlappingExomeIntervals.size())});

        return result.toArray(new Object[result.size()][]);
    }

    @DataProvider(name="correctRangeObjectInitializationData")
    public Object[][] correctRangeObjectInitializationData() {
        return new Object[][] {
                { 0, 0},
                { 0, 1},
                { 0, 100},
                { 10, 10}
        };
    }

    @DataProvider(name="wrongRangeObjectInitializationData")
    public Object[][] wrongRangeObjectInitializationData() {
        return new Object[][] {
                { -1 , 0},
                { 20, 19},
        };
    }

    @DataProvider(name="invalidRangeData")
    public Object[][] invalidRangeData() {
        return new Object[][] {
                { -1 , 0},
                { 2, 0},
                { nonOverlappingExomeIntervals.size() + 1, nonOverlappingExomeIntervals.size() }
        };
    }

    @DataProvider(name="arbitraryExonRageData")
    public Object[][] arbitraryExonRageData() {
        final int exonCount = nonOverlappingExomeIntervals.size();
        final Random rdn = new Random(1313); // feeling very lucky.
        final Object[][] result = new Object[100][];
        // make sure there is at least a one element exon range request.
        result[0] = new Object[] { exonCount - 4, exonCount - 3 };
        for (int i = 1; i < result.length; i++) {
            final int from = rdn.nextInt(exonCount);
            final int to = rdn.nextInt(exonCount - from) + from;
            result[i] = new Object[] { from , to};
        }
        return result;
    }

    @Test
    public void testRandomLookupExonByLocation() {
        final Object[][] lookups = exonLookUpData();
        final Random rdn = new Random(11132331);
        final IntervalBackedExonCollection exonDB = new IntervalBackedExonCollection(nonOverlappingExomeIntervals);
        for (int i = 0; i < 1000; i++) {
            final Object[] params = lookups[rdn.nextInt(lookups.length)];
            @SuppressWarnings("unchecked")
            final GenomeLoc query = (GenomeLoc) params[0];
            @SuppressWarnings("unchecked")
            final GenomeLoc expected = (GenomeLoc) params[1];
            final GenomeLoc observed = exonDB.exon(query);
            Assert.assertEquals(observed,expected);
        }

    }

    @DataProvider(name="exonLookUpData")
    public Object[][] exonLookUpData() {
        final List<Object[]> result = new ArrayList<>();
        for (int i = 0; i < nonOverlappingExomeIntervals.size(); i++) {
            result.add(new Object[]{nonOverlappingExomeIntervals.get(i), nonOverlappingExomeIntervals.get(i), i});
        }
        for (int i = 0; i < nonOverlappingExomeIntervals.size(); i++) {
            result.add(new Object[]{nonOverlappingExomeIntervals.get(i).getStartLocation(), nonOverlappingExomeIntervals.get(i), i});
        }
        for (int i = 0; i < nonOverlappingExomeIntervals.size(); i++) {
            result.add(new Object[]{nonOverlappingExomeIntervals.get(i).getStopLocation(), nonOverlappingExomeIntervals.get(i), i});
        }
        for (int i = 1; i < nonOverlappingExomeIntervals.size(); i++) {
            final GenomeLoc previous = nonOverlappingExomeIntervals.get(i-1);
            final GenomeLoc next = nonOverlappingExomeIntervals.get(i);
            final GenomeLoc query = previous.getContig().equals(next.getContig())
                    ? ExomeToolsTestUtils.createGenomeLoc(previous.getContig(),previous.getStop() + 1,next.getStart() - 1)
                    : ExomeToolsTestUtils.createGenomeLoc(next.getContig(),1,next.getStart() - 1);
            result.add(new Object[] { query , null, -i-1});
        }
        return result.toArray(new Object[result.size()][]);
    }

    @DataProvider(name="wrongExonLookUpData")
    public Object[][] wrongExonLookUpData() {
        return new Object[][] {
                {null},
                {ExomeToolsTestUtils.createOverEntireContig(nonOverlappingExomeIntervals.get(0).getContig())},
                {ExomeToolsTestUtils.createGenomeLoc(nonOverlappingExomeIntervals.get(0).getContig(),
                        nonOverlappingExomeIntervals.get(0).getStop(),
                        nonOverlappingExomeIntervals.get(1).getStart())}
        };
    }
}
