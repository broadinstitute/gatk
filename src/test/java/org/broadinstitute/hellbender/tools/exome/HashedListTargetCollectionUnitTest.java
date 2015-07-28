package org.broadinstitute.hellbender.tools.exome;

import htsjdk.samtools.SAMSequenceRecord;
import org.broadinstitute.hellbender.utils.IndexRange;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.test.ExomeToolsTestUtils;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;
import java.util.stream.Collectors;

/**
 * Unit tests for {@link HashedListTargetCollection}.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class HashedListTargetCollectionUnitTest extends BaseTest {

    /**
     * Average target size for randomly generated test data.
     */
    private static final int averageExonSize = 100;

    /**
     * Std. dev. for target size for randomly generated test data.
     */
    private static final float sdExonSize = 20;

    /**
     * Average inter-target gap in bp for randomly generated test data.
     */
    private static final int averageExonIntergapSize = 100;

    /**
     * Std. dev. for inter-target gap in bp for randomly generated test data.
     */
    private static final float sdExonIntergapSize = 20;

    /**
     * Minimum inter-target gap size for randomly generated test data.
     */
    private static final int minimumExonIntergapSize = 20;

    /**
     * Minimum target size for randomly generated test data.
     */
    private static final int minimumExonSize = 1;

    /**
     * Reference to the test intervals.
     * <p>
     *     Initialized before all tests in {@link #setUp}.
     * </p>
     */
    private List<SimpleInterval> nonOverlappingExomeIntervals;

    /**
     * A HashedListTargetCollection of the test intervals.
     * <p>
     *     Initialized before all tests in {@link #setUp}.
     * </p>
     */
    private HashedListTargetCollection<SimpleInterval> exonDB;

    @BeforeClass
    public void setUp() {
        nonOverlappingExomeIntervals = new ArrayList<>(100);
        final Random rdn = new Random(13);// some "random" but fixed seed to make sure errors are deterministic.
        for (int i = 0; i < ExomeToolsTestUtils.REFERENCE_DICTIONARY.size(); i++) {
            int current = 0;
            final SAMSequenceRecord sequence = ExomeToolsTestUtils.REFERENCE_DICTIONARY.getSequence(i);
            while (current < sequence.getSequenceLength()) {
                int start = current + Math.max(minimumExonIntergapSize,
                        (int) Math.round(rdn.nextGaussian() * sdExonIntergapSize + averageExonIntergapSize));
                if (start >= sequence.getSequenceLength()) {
                    break;
                }
                int size = Math.max(minimumExonSize,
                        (int) Math.round(rdn.nextGaussian() * sdExonSize + averageExonSize));
                int stop = start + size - 1;
                if (stop >= sequence.getSequenceLength()) {
                    break;
                }
                nonOverlappingExomeIntervals.add(
                        ExomeToolsTestUtils.createInterval(
                                sequence.getSequenceName(), start, stop));
                current = stop + 1;
            }
        }
        Collections.sort(nonOverlappingExomeIntervals, IntervalUtils.LEXICOGRAPHICAL_ORDER_COMPARATOR);
        exonDB = new HashedListTargetCollection<>(nonOverlappingExomeIntervals);
    }

    @Test
    public void testCorrectInitialization() {
        new HashedListTargetCollection<>(nonOverlappingExomeIntervals);
    }

    @Test
    public void testCorrectInitializationSortingContigsByName() {
        final List<SimpleInterval> sortedLexicographically = nonOverlappingExomeIntervals.stream()
                .sorted(IntervalUtils.LEXICOGRAPHICAL_ORDER_COMPARATOR).collect(Collectors.toList());
        final HashedListTargetCollection<SimpleInterval> exonCollection = new HashedListTargetCollection<>(nonOverlappingExomeIntervals);
        Assert.assertEquals(exonCollection.targetCount(),nonOverlappingExomeIntervals.size());
        for (int i = 0; i < nonOverlappingExomeIntervals.size(); i++) {
            Assert.assertEquals(exonCollection.target(i),sortedLexicographically.get(i));
        }
    }

    @Test
    public void testCorrectInitializationUnsorted() {
        final List<SimpleInterval> unsortedIntervals = new ArrayList<>(nonOverlappingExomeIntervals);
        final Random rdn = new Random(13313);
        // shuffle forward:
        for (int i = 0; i < unsortedIntervals.size() - 1; i++) {
            final int otherIndex = rdn.nextInt(unsortedIntervals.size() - i - 1) + i + 1;
            // do a swap of elements i <-> otherIndex:
            unsortedIntervals.set(i,unsortedIntervals.set(otherIndex,unsortedIntervals.get(i)));
        }
        // shuffle backwards:
        for (int i = unsortedIntervals.size() - 1; i > 0; --i) {
            final int otherIndex = rdn.nextInt(i);
            // do a swap of elements i <-> otherIndex:
            unsortedIntervals.set(i,unsortedIntervals.set(otherIndex,unsortedIntervals.get(i)));
        }
        // make sure that there are at least two intervals out of order.

        boolean foundOutOfOrder = false;
        OUTER_LOOP: for (int i = 0; i < unsortedIntervals.size(); i++) {
            final SimpleInterval iInterval = unsortedIntervals.get(i);
            for (int j = i + 1; j < unsortedIntervals.size(); j++) {
                final SimpleInterval jInterval = unsortedIntervals.get(j);
                if (IntervalUtils.LEXICOGRAPHICAL_ORDER_COMPARATOR.compare(iInterval, jInterval) > 0) {
                    foundOutOfOrder = true;
                    break OUTER_LOOP;
                }
            }
        }
        Assert.assertTrue(foundOutOfOrder,"the order randomization step is not working");
        final HashedListTargetCollection<SimpleInterval> exonCollection = new HashedListTargetCollection<>(unsortedIntervals);
        Assert.assertEquals(exonCollection.targetCount(),nonOverlappingExomeIntervals.size());
        for (int i = 0; i < nonOverlappingExomeIntervals.size(); i++) {
            Assert.assertEquals(exonCollection.target(i),nonOverlappingExomeIntervals.get(i));
        }
    }

    @Test(dataProvider="correctRangeObjectInitializationData")
    public void testCorrectRangeObjectInitialization(final int from, final int to) {
        new IndexRange(from,to);
    }

    @Test
    public void testCorrectEmptyInitialization() {
        new HashedListTargetCollection<>(Collections.emptyList());
    }

    @Test(expectedExceptions = {IllegalArgumentException.class})
    public void testIncorrectInitialization1() {
        new HashedListTargetCollection<>(null);
    }

    @Test(expectedExceptions = {IllegalArgumentException.class})
    public void testIncorrectInitialization2() {
        final List<SimpleInterval> intervals = new ArrayList<>(nonOverlappingExomeIntervals);
        intervals.add(nonOverlappingExomeIntervals.size() >> 1,null);
        new HashedListTargetCollection<>(intervals);
    }

    @Test(expectedExceptions = {IllegalArgumentException.class})
    public void testIncorrectInitialization5() {
        final List<SimpleInterval> intervals = new ArrayList<>(nonOverlappingExomeIntervals);
        intervals.add(1,intervals.get(0));
        new HashedListTargetCollection<>(intervals);
    }

    @Test(expectedExceptions = {IllegalArgumentException.class})
    public void testIncorrectInitialization6() {
        final List<SimpleInterval> intervals = new ArrayList<>(nonOverlappingExomeIntervals);
        final SimpleInterval interval = intervals.get(0);
        intervals.add(1, new SimpleInterval(interval));
        new HashedListTargetCollection<>(intervals);
    }

    @Test(dependsOnMethods={"testCorrectInitialization"})
    public void testExonCount() {
        Assert.assertEquals(exonDB.targetCount(),nonOverlappingExomeIntervals.size());
    }

    @Test(dependsOnMethods={"testCorrectInitialization"})
    public void testExomeSize() {
        int size = 0;
        for (final SimpleInterval loc : nonOverlappingExomeIntervals) {
            size += loc.size();
        }
        Assert.assertEquals(exonDB.exomeSize(),size);
    }

    @Test(dependsOnMethods = {"testCorrectInitialization"})
    public void testExons() {
        final List<SimpleInterval> exons = exonDB.targets();
        Assert.assertNotNull(exons);
        Assert.assertEquals(exons.size(),nonOverlappingExomeIntervals.size());
        for (int i = 0; i < exons.size(); i++) {
            Assert.assertEquals(exons.get(i), nonOverlappingExomeIntervals.get(i));
        }
    }

    @Test(dependsOnMethods = {"testCorrectInitialization"})
    public void testExonsFullRange() {
        final List<SimpleInterval> exons = exonDB.targets(0, exonDB.targetCount());
        Assert.assertNotNull(exons);
        Assert.assertEquals(exons.size(),nonOverlappingExomeIntervals.size());
        for (int i = 0; i < exons.size(); i++) {
            Assert.assertEquals(exons.get(i), nonOverlappingExomeIntervals.get(i));
        }
    }

    @Test(dependsOnMethods = {"testCorrectInitialization","testCorrectRangeObjectInitialization"})
    public void testExonFullRangeObject() {
        final IndexRange range = new IndexRange(0,exonDB.targetCount());
        final List<SimpleInterval> exons = exonDB.targets(range);
        Assert.assertNotNull(exons);
        Assert.assertEquals(exons.size(),nonOverlappingExomeIntervals.size());
        for (int i = 0; i < exons.size(); i++) {
            Assert.assertEquals(exons.get(i), nonOverlappingExomeIntervals.get(i));
        }
    }

    @Test(dependsOnMethods = {"testCorrectInitialization"})
    public void testEmptyRange() {
        final List<SimpleInterval> exons = exonDB.targets(0, 0);
        Assert.assertNotNull(exons);
        Assert.assertEquals(exons.size(),0);
    }

    @Test(dependsOnMethods = {"testCorrectInitialization","testCorrectRangeObjectInitialization"})
    public void testEmptyRangeObject() {
        final IndexRange range = new IndexRange(0,0);
        final List<SimpleInterval> exons = exonDB.targets(range);
        Assert.assertNotNull(exons);
        Assert.assertEquals(exons.size(),0);
    }

    @Test(dependsOnMethods = {"testCorrectInitialization"},
          dataProvider = "arbitraryExonRageData")
    public void testArbitraryExonRange(final int from, final int to) {
        final List<SimpleInterval> exons = exonDB.targets(from, to);
        Assert.assertNotNull(exons);
        final int expectedSize = to - from;
        Assert.assertEquals(exons.size(), expectedSize);
        for (int i = from; i < to; i++) {
            Assert.assertEquals(exons.get(i - from), nonOverlappingExomeIntervals.get(i));
        }
    }

    @Test(dependsOnMethods = {"testCorrectInitialization","testCorrectRangeObjectInitialization"},
            dataProvider = "arbitraryExonRageData")
    public void testArbitraryExonRangeObject(final int from, final int to) {
        final IndexRange range = new IndexRange(from,to);
        final List<SimpleInterval> exons = exonDB.targets(range);
        Assert.assertNotNull(exons);
        final int expectedSize = to - from;
        Assert.assertEquals(exons.size(), expectedSize);
        for (int i = from; i < to; i++) {
            Assert.assertEquals(exons.get(i - from), nonOverlappingExomeIntervals.get(i));
        }
    }

    @Test(dependsOnMethods = "testCorrectInitialization")
    public void testLocation() {
        for (final SimpleInterval loc : nonOverlappingExomeIntervals) {
            Assert.assertEquals(exonDB.location(loc), loc);
        }
    }

    @Test(dependsOnMethods = "testCorrectInitialization",
        expectedExceptions = IllegalArgumentException.class )
    public void testWrongLocation() {
        exonDB.location(null);
    }

    @Test(dependsOnMethods = "testCorrectInitialization")
    public void testLocationFromIndex() {
        for (int i = 0; i < nonOverlappingExomeIntervals.size(); i++) {
            Assert.assertEquals(exonDB.location(i), nonOverlappingExomeIntervals.get(i));
        }
    }

    @Test(dependsOnMethods = {"testCorrectInitialization"},
          expectedExceptions = {IndexOutOfBoundsException.class, IllegalArgumentException.class},
          dataProvider = "invalidRangeData")
    public void testInvalidRange(final int from, final int to) {
        exonDB.targets(from, to);
    }

    @Test(dependsOnMethods = {"testCorrectInitialization","testCorrectRangeObjectInitialization"},
            expectedExceptions = {IndexOutOfBoundsException.class, IllegalArgumentException.class},
            dataProvider = "invalidRangeData")
    public void testInvalidRangeUsingObject(final int from, final int to) {
        final IndexRange range = new IndexRange(from,to);
        exonDB.targets(range);
    }

    @Test(dependsOnMethods = "testCorrectInitialization",
          dataProvider = "exonLookUpData")
    public void testExonByLocation(final SimpleInterval location, final SimpleInterval expected,
                                   @SuppressWarnings("unused") final int index) {
        final SimpleInterval actual = exonDB.target(location);
        if (expected == null) {
            Assert.assertNull(actual);
        } else {
            Assert.assertNotNull(actual);
            Assert.assertEquals(actual,expected);
        }
    }

    @Test(dependsOnMethods = "testCorrectInitialization",
              dataProvider = "exonLookUpData")
    public void testExonIndexByLocation(final SimpleInterval location, @SuppressWarnings("unused")  final SimpleInterval expected,
                                   final int index) {
        final int actual = exonDB.index(location);
        Assert.assertEquals(actual, index);
    }

    @Test(dependsOnMethods = "testCorrectInitialization",
              dataProvider = "wrongExonLookUpData",
            expectedExceptions = {TargetCollection.AmbiguousTargetException.class, IllegalArgumentException.class})
    public void testWrongExonIndexByLocation(final SimpleInterval location) {
        exonDB.index(location);
    }

    @Test(dependsOnMethods = "testCorrectInitialization",
            dataProvider = "wrongExonLookUpData",
            expectedExceptions = {TargetCollection.AmbiguousTargetException.class, IllegalArgumentException.class})
    public void testWrongExonByLocation(final SimpleInterval location) {
        exonDB.target(location);
    }

    @Test(dependsOnMethods = "testCorrectInitialization")
    public void testExonByName() {
        for (final SimpleInterval loc : nonOverlappingExomeIntervals) {
            Assert.assertEquals(exonDB.target(loc.toString()), loc);
        }
        Assert.assertNull(exonDB.target("no-id"));
    }

    @Test(dependsOnMethods = "testCorrectInitialization")
    public void testExonIndexByName() {
        for (int i = 0; i < nonOverlappingExomeIntervals.size(); i++) {
            Assert.assertEquals(exonDB.index(nonOverlappingExomeIntervals.get(i)),i);
        }
        Assert.assertEquals(exonDB.index("no-id"), -1);
    }

    @Test(dependsOnMethods = {"testCorrectInitialization"},
        dataProvider="indexRangeData")
    public void testIndexRange(final SimpleInterval region, final IndexRange expected) {
        final IndexRange range = exonDB.indexRange(region);
        Assert.assertEquals(range, expected);
        Assert.assertEquals(exonDB.targetCount(region), range.size());
    }

    @Test(dependsOnMethods = {"testCorrectInitialization"},
            dataProvider="indexRangeData")
    public void testForEachExon(final SimpleInterval region, final IndexRange expected) {
        final List<SimpleInterval> exons = exonDB.targets(region);
        final List<Integer> exonIndices = new ArrayList<>(exons.size());
        final List<SimpleInterval> exonObjects = new ArrayList<>(exons.size());
        exonDB.forEachTarget(region, (index, exon) -> {
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
    public void testExons(final SimpleInterval region, final IndexRange expected) {
        final List<SimpleInterval> exons = exonDB.targets(region);
        final List<SimpleInterval> expectedExons =  nonOverlappingExomeIntervals.subList(expected.from,expected.to);
        Assert.assertEquals(exons,expectedExons);
    }

    @Test(dependsOnMethods={"testCorrectInitialization","testExonCount"})
    public void testExonName() {
        for (int i = 0; i < exonDB.targetCount(); i++) {
            final SimpleInterval loc = exonDB.target(i);
            final String name = exonDB.name(loc);
            final SimpleInterval observed = exonDB.target(name);
            Assert.assertEquals(observed,loc);
            Assert.assertEquals(exonDB.index(name),i);
        }
        Assert.assertNull(exonDB.target("no-interval"));
        Assert.assertEquals(exonDB.index("no-interval"),-1);
    }

    @DataProvider(name="indexRangeData")
    public Object[][] indexRangeData() {
        final List<Object[]> result = new ArrayList<>();

        final Random rdn = new Random(131313);
        for (int i = 0; i < nonOverlappingExomeIntervals.size(); i++) {
            result.add(new Object[]{nonOverlappingExomeIntervals.get(i), new IndexRange(i, i + 1)});
            result.add(new Object[]{
                    ExomeToolsTestUtils.createInterval(nonOverlappingExomeIntervals.get(i).getContig(),
                            nonOverlappingExomeIntervals.get(i).getStart() - 1), new IndexRange(i,i)});
            int j;
            for (j = i + 1; j < nonOverlappingExomeIntervals.size() && rdn.nextBoolean(); j++) {
                if (!nonOverlappingExomeIntervals.get(j).getContig().equals(nonOverlappingExomeIntervals.get(i).getContig())) {
                    break;
                }
            }
            result.add(new Object[]{
                    ExomeToolsTestUtils.createInterval(nonOverlappingExomeIntervals.get(i).getContig(),
                            nonOverlappingExomeIntervals.get(i).getStart(),
                            nonOverlappingExomeIntervals.get(j - 1).getEnd()),
                    new IndexRange(i, j)
            });
        }
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
        // make sure there is at least a one element target range request.
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
        for (int i = 0; i < 1000; i++) {
            final Object[] params = lookups[rdn.nextInt(lookups.length)];
            @SuppressWarnings("unchecked")
            final SimpleInterval query = (SimpleInterval) params[0];
            @SuppressWarnings("unchecked")
            final SimpleInterval expected = (SimpleInterval) params[1];
            final SimpleInterval observed = exonDB.target(query);
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
            final SimpleInterval interval = nonOverlappingExomeIntervals.get(i);
            result.add(new Object[]{new SimpleInterval(interval.getContig(),interval.getStart(),interval.getStart()), interval, i});
        }
        for (int i = 0; i < nonOverlappingExomeIntervals.size(); i++) {
            final SimpleInterval interval = nonOverlappingExomeIntervals.get(i);
            result.add(new Object[]{new SimpleInterval(interval.getContig(),interval.getEnd(),interval.getEnd()), interval, i});
        }
        for (int i = 1; i < nonOverlappingExomeIntervals.size(); i++) {
            final SimpleInterval previous = nonOverlappingExomeIntervals.get(i-1);
            final SimpleInterval next = nonOverlappingExomeIntervals.get(i);
            final SimpleInterval query = previous.getContig().equals(next.getContig())
                    ? ExomeToolsTestUtils.createInterval(previous.getContig(), previous.getEnd() + 1, next.getStart() - 1)
                    : ExomeToolsTestUtils.createInterval(next.getContig(), 1, next.getStart() - 1);
            result.add(new Object[]{query, null, -i - 1});
        }
        return result.toArray(new Object[result.size()][]);
    }

    @DataProvider(name="wrongExonLookUpData")
    public Object[][] wrongExonLookUpData() {
        return new Object[][] {
                {null},
                {ExomeToolsTestUtils.createOverEntireContig(nonOverlappingExomeIntervals.get(0).getContig())},
                {ExomeToolsTestUtils.createInterval(nonOverlappingExomeIntervals.get(0).getContig(),
                        nonOverlappingExomeIntervals.get(0).getEnd(),
                        nonOverlappingExomeIntervals.get(1).getStart())}
        };
    }

    @DataProvider(name="trimmedIntervalData")
    public Object[][] trimmedIntervalData() {

        final List<SimpleInterval> intervals = nonOverlappingExomeIntervals;    //alias to reduce clutter

        final String contig = intervals.get(0).getContig();
        //find the last index of the first contig
        int j = 0;
        while (j < intervals.size() && intervals.get(j).getContig().equals(contig)) {
            j++;
        }
        j--;

        //shift of untrimmed interval ends
        final int offset = minimumExonIntergapSize/2;

         // Each data row is {untrimmed interval, epected trimmed interval
        final int firstTargetStart = intervals.get(0).getStart();
        final int firstTargetEnd = intervals.get(0).getEnd();
        final int lastTargetEnd = intervals.get(j).getEnd();
        return new Object[][] {
                //bracket a single target
                {new SimpleInterval(contig, firstTargetStart - offset, firstTargetEnd + offset), new SimpleInterval(contig, firstTargetStart, firstTargetEnd)},
                //within a single target
                {new SimpleInterval(contig, firstTargetEnd -1 , firstTargetEnd + offset), new SimpleInterval(contig, firstTargetStart, firstTargetEnd)},
                //a lot of targets
                {new SimpleInterval(contig, firstTargetStart - offset, lastTargetEnd + offset), new SimpleInterval(contig, firstTargetStart, lastTargetEnd)},
                //past the last target
                {new SimpleInterval(contig, lastTargetEnd + 1, lastTargetEnd + 2), new SimpleInterval(contig, lastTargetEnd + 1, lastTargetEnd + 1)},
                //a non-existant contig
                {new SimpleInterval("non-existant contig", 1, 10), new SimpleInterval("non-existant contig", 1, 1)}
        };
    }

    @Test(dataProvider = "trimmedIntervalData")
    public void testTrimmedInterval(final SimpleInterval untrimmed, final SimpleInterval trimmed) {
        Assert.assertEquals(exonDB.trimmedInterval(untrimmed), trimmed);
    }
}
