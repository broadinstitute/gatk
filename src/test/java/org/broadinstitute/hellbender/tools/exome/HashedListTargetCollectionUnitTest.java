package org.broadinstitute.hellbender.tools.exome;

import htsjdk.samtools.SAMSequenceRecord;
import org.broadinstitute.hellbender.utils.IndexRange;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.test.TargetsToolsTestUtils;
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
    private static final int averageTargetSize = 100;

    /**
     * Std. dev. for target size for randomly generated test data.
     */
    private static final float sdTargetSize = 20;

    /**
     * Average inter-target gap in bp for randomly generated test data.
     */
    private static final int averageTargetIntergapSize = 100;

    /**
     * Std. dev. for inter-target gap in bp for randomly generated test data.
     */
    private static final float sdTargetIntergapSize = 20;

    /**
     * Minimum inter-target gap size for randomly generated test data.
     */
    private static final int minimumTargetIntergapSize = 20;

    /**
     * Minimum target size for randomly generated test data.
     */
    private static final int minimumTargetSize = 1;

    /**
     * Reference to the test intervals.
     * <p>
     *     Initialized before all tests in {@link #setUp}.
     * </p>
     */
    private List<SimpleInterval> nonOverlappingTargetIntervals;

    /**
     * A HashedListTargetCollection of the test intervals.
     * <p>
     *     Initialized before all tests in {@link #setUp}.
     * </p>
     */
    private HashedListTargetCollection<SimpleInterval> targetDB;

    @BeforeClass
    public void setUp() {
        nonOverlappingTargetIntervals = new ArrayList<>(100);
        final Random rdn = new Random(13);// some "random" but fixed seed to make sure errors are deterministic.
        for (int i = 0; i < TargetsToolsTestUtils.REFERENCE_DICTIONARY.size(); i++) {
            int current = 0;
            final SAMSequenceRecord sequence = TargetsToolsTestUtils.REFERENCE_DICTIONARY.getSequence(i);
            while (current < sequence.getSequenceLength()) {
                int start = current + Math.max(minimumTargetIntergapSize,
                        (int) Math.round(rdn.nextGaussian() * sdTargetIntergapSize + averageTargetIntergapSize));
                if (start >= sequence.getSequenceLength()) {
                    break;
                }
                int size = Math.max(minimumTargetSize,
                        (int) Math.round(rdn.nextGaussian() * sdTargetSize + averageTargetSize));
                int stop = start + size - 1;
                if (stop >= sequence.getSequenceLength()) {
                    break;
                }
                nonOverlappingTargetIntervals.add(
                        TargetsToolsTestUtils.createInterval(
                                sequence.getSequenceName(), start, stop));
                current = stop + 1;
            }
        }
        Collections.sort(nonOverlappingTargetIntervals, IntervalUtils.LEXICOGRAPHICAL_ORDER_COMPARATOR);
        targetDB = new HashedListTargetCollection<>(nonOverlappingTargetIntervals);
    }

    @Test
    public void testCorrectInitialization() {
        new HashedListTargetCollection<>(nonOverlappingTargetIntervals);
    }

    @Test
    public void testCorrectInitializationSortingContigsByName() {
        final List<SimpleInterval> sortedLexicographically = nonOverlappingTargetIntervals.stream()
                .sorted(IntervalUtils.LEXICOGRAPHICAL_ORDER_COMPARATOR).collect(Collectors.toList());
        final HashedListTargetCollection<SimpleInterval> targetCollection = new HashedListTargetCollection<>(nonOverlappingTargetIntervals);
        Assert.assertEquals(targetCollection.targetCount(), nonOverlappingTargetIntervals.size());
        for (int i = 0; i < nonOverlappingTargetIntervals.size(); i++) {
            Assert.assertEquals(targetCollection.target(i),sortedLexicographically.get(i));
        }
    }

    @Test
    public void testCorrectInitializationUnsorted() {
        final List<SimpleInterval> unsortedIntervals = new ArrayList<>(nonOverlappingTargetIntervals);
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
        final HashedListTargetCollection<SimpleInterval> targetCollection = new HashedListTargetCollection<>(unsortedIntervals);
        Assert.assertEquals(targetCollection.targetCount(), nonOverlappingTargetIntervals.size());
        for (int i = 0; i < nonOverlappingTargetIntervals.size(); i++) {
            Assert.assertEquals(targetCollection.target(i), nonOverlappingTargetIntervals.get(i));
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
        final List<SimpleInterval> intervals = new ArrayList<>(nonOverlappingTargetIntervals);
        intervals.add(nonOverlappingTargetIntervals.size() >> 1,null);
        new HashedListTargetCollection<>(intervals);
    }

    @Test(expectedExceptions = {IllegalArgumentException.class})
    public void testIncorrectInitialization5() {
        final List<SimpleInterval> intervals = new ArrayList<>(nonOverlappingTargetIntervals);
        intervals.add(1,intervals.get(0));
        new HashedListTargetCollection<>(intervals);
    }

    @Test(expectedExceptions = {IllegalArgumentException.class})
    public void testIncorrectInitialization6() {
        final List<SimpleInterval> intervals = new ArrayList<>(nonOverlappingTargetIntervals);
        final SimpleInterval interval = intervals.get(0);
        intervals.add(1, new SimpleInterval(interval));
        new HashedListTargetCollection<>(intervals);
    }

    @Test(dependsOnMethods={"testCorrectInitialization"})
    public void testTargetCount() {
        Assert.assertEquals(targetDB.targetCount(), nonOverlappingTargetIntervals.size());
    }

    @Test(dependsOnMethods={"testCorrectInitialization"})
    public void testTotalSize() {
        int size = 0;
        for (final SimpleInterval loc : nonOverlappingTargetIntervals) {
            size += loc.size();
        }
        Assert.assertEquals(targetDB.totalSize(),size);
    }

    @Test(dependsOnMethods = {"testCorrectInitialization"})
    public void testTargets() {
        final List<SimpleInterval> targets = targetDB.targets();
        Assert.assertNotNull(targets);
        Assert.assertEquals(targets.size(), nonOverlappingTargetIntervals.size());
        for (int i = 0; i < targets.size(); i++) {
            Assert.assertEquals(targets.get(i), nonOverlappingTargetIntervals.get(i));
        }
    }

    @Test(dependsOnMethods = {"testCorrectInitialization"})
    public void testTargetsFullRange() {
        final List<SimpleInterval> targets = targetDB.targets(0, targetDB.targetCount());
        Assert.assertNotNull(targets);
        Assert.assertEquals(targets.size(), nonOverlappingTargetIntervals.size());
        for (int i = 0; i < targets.size(); i++) {
            Assert.assertEquals(targets.get(i), nonOverlappingTargetIntervals.get(i));
        }
    }

    @Test(dependsOnMethods = {"testCorrectInitialization","testCorrectRangeObjectInitialization"})
    public void testTargetFullRangeObject() {
        final IndexRange range = new IndexRange(0, targetDB.targetCount());
        final List<SimpleInterval> targets = targetDB.targets(range);
        Assert.assertNotNull(targets);
        Assert.assertEquals(targets.size(), nonOverlappingTargetIntervals.size());
        for (int i = 0; i < targets.size(); i++) {
            Assert.assertEquals(targets.get(i), nonOverlappingTargetIntervals.get(i));
        }
    }

    @Test(dependsOnMethods = {"testCorrectInitialization"})
    public void testEmptyRange() {
        final List<SimpleInterval> targets = targetDB.targets(0, 0);
        Assert.assertNotNull(targets);
        Assert.assertEquals(targets.size(),0);
    }

    @Test(dependsOnMethods = {"testCorrectInitialization","testCorrectRangeObjectInitialization"})
    public void testEmptyRangeObject() {
        final IndexRange range = new IndexRange(0,0);
        final List<SimpleInterval> targets = targetDB.targets(range);
        Assert.assertNotNull(targets);
        Assert.assertEquals(targets.size(),0);
    }

    @Test(dependsOnMethods = {"testCorrectInitialization"},
          dataProvider = "arbitraryTargetRageData")
    public void testArbitraryTargetRange(final int from, final int to) {
        final List<SimpleInterval> targets = targetDB.targets(from, to);
        Assert.assertNotNull(targets);
        final int expectedSize = to - from;
        Assert.assertEquals(targets.size(), expectedSize);
        for (int i = from; i < to; i++) {
            Assert.assertEquals(targets.get(i - from), nonOverlappingTargetIntervals.get(i));
        }
    }

    @Test(dependsOnMethods = {"testCorrectInitialization","testCorrectRangeObjectInitialization"},
            dataProvider = "arbitraryTargetRageData")
    public void testArbitraryTargetRangeObject(final int from, final int to) {
        final IndexRange range = new IndexRange(from,to);
        final List<SimpleInterval> targets = targetDB.targets(range);
        Assert.assertNotNull(targets);
        final int expectedSize = to - from;
        Assert.assertEquals(targets.size(), expectedSize);
        for (int i = from; i < to; i++) {
            Assert.assertEquals(targets.get(i - from), nonOverlappingTargetIntervals.get(i));
        }
    }

    @Test(dependsOnMethods = "testCorrectInitialization")
    public void testLocation() {
        for (final SimpleInterval loc : nonOverlappingTargetIntervals) {
            Assert.assertEquals(targetDB.location(loc), loc);
        }
    }

    @Test(dependsOnMethods = "testCorrectInitialization",
        expectedExceptions = IllegalArgumentException.class )
    public void testWrongLocation() {
        targetDB.location(null);
    }

    @Test(dependsOnMethods = "testCorrectInitialization")
    public void testLocationFromIndex() {
        for (int i = 0; i < nonOverlappingTargetIntervals.size(); i++) {
            Assert.assertEquals(targetDB.location(i), nonOverlappingTargetIntervals.get(i));
        }
    }

    @Test(dependsOnMethods = {"testCorrectInitialization"},
          expectedExceptions = {IndexOutOfBoundsException.class, IllegalArgumentException.class},
          dataProvider = "invalidRangeData")
    public void testInvalidRange(final int from, final int to) {
        targetDB.targets(from, to);
    }

    @Test(dependsOnMethods = {"testCorrectInitialization","testCorrectRangeObjectInitialization"},
            expectedExceptions = {IndexOutOfBoundsException.class, IllegalArgumentException.class},
            dataProvider = "invalidRangeData")
    public void testInvalidRangeUsingObject(final int from, final int to) {
        final IndexRange range = new IndexRange(from,to);
        targetDB.targets(range);
    }

    @Test(dependsOnMethods = "testCorrectInitialization",
          dataProvider = "targetLookUpData")
    public void testTargetByLocation(final SimpleInterval location, final SimpleInterval expected,
                                     @SuppressWarnings("unused") final int index) {
        final SimpleInterval actual = targetDB.target(location);
        if (expected == null) {
            Assert.assertNull(actual);
        } else {
            Assert.assertNotNull(actual);
            Assert.assertEquals(actual,expected);
        }
    }

    @Test(dependsOnMethods = "testCorrectInitialization",
              dataProvider = "targetLookUpData")
    public void testTargetIndexByLocation(final SimpleInterval location, @SuppressWarnings("unused") final SimpleInterval expected,
                                          final int index) {
        final int actual = targetDB.index(location);
        Assert.assertEquals(actual, index);
    }

    @Test(dependsOnMethods = "testCorrectInitialization",
              dataProvider = "wrongTargetLookUpData",
            expectedExceptions = {TargetCollection.AmbiguousTargetException.class, IllegalArgumentException.class})
    public void testWrongTargetIndexByLocation(final SimpleInterval location) {
        targetDB.index(location);
    }

    @Test(dependsOnMethods = "testCorrectInitialization",
            dataProvider = "wrongTargetLookUpData",
            expectedExceptions = {TargetCollection.AmbiguousTargetException.class, IllegalArgumentException.class})
    public void testWrongTargetByLocation(final SimpleInterval location) {
        targetDB.target(location);
    }

    @Test(dependsOnMethods = "testCorrectInitialization")
    public void testTargetByName() {
        for (final SimpleInterval loc : nonOverlappingTargetIntervals) {
            Assert.assertEquals(targetDB.target(loc.toString()), loc);
        }
        Assert.assertNull(targetDB.target("no-id"));
    }

    @Test(dependsOnMethods = "testCorrectInitialization")
    public void testTargetIndexByName() {
        for (int i = 0; i < nonOverlappingTargetIntervals.size(); i++) {
            Assert.assertEquals(targetDB.index(nonOverlappingTargetIntervals.get(i)),i);
        }
        Assert.assertEquals(targetDB.index("no-id"), -1);
    }

    @Test(dependsOnMethods = {"testCorrectInitialization"},
        dataProvider="indexRangeData")
    public void testIndexRange(final SimpleInterval region, final IndexRange expected) {
        final IndexRange range = targetDB.indexRange(region);
        Assert.assertEquals(range, expected);
        Assert.assertEquals(targetDB.targetCount(region), range.size());
    }

    @Test(dependsOnMethods = {"testCorrectInitialization"},
            dataProvider="indexRangeData")
    public void testForEachTarget(final SimpleInterval region, final IndexRange expected) {
        final List<SimpleInterval> targets = targetDB.targets(region);
        final List<Integer> targetIndices = new ArrayList<>(targets.size());
        final List<SimpleInterval> targetObjects = new ArrayList<>(targets.size());
        targetDB.forEachTarget(region, (index, target) -> {
            targetIndices.add(index);
            targetObjects.add(target);
        });
        Assert.assertEquals(targetIndices.size(),targets.size());
        for (int i = 0; i < targets.size(); i++) {
            Assert.assertEquals(targetIndices.get(i).intValue(),expected.from + i);
            Assert.assertEquals(targetObjects.get(i),targets.get(i));
        }
    }

    @Test(dependsOnMethods = {"testCorrectInitialization"},
       dataProvider = "indexRangeData")
    public void testTargets(final SimpleInterval region, final IndexRange expected) {
        final List<SimpleInterval> targets = targetDB.targets(region);
        final List<SimpleInterval> expectedTargets =  nonOverlappingTargetIntervals.subList(expected.from,expected.to);
        Assert.assertEquals(targets,expectedTargets);
    }

    @Test(dependsOnMethods={"testCorrectInitialization", "testTargetCount"})
    public void testTargetName() {
        for (int i = 0; i < targetDB.targetCount(); i++) {
            final SimpleInterval loc = targetDB.target(i);
            final String name = targetDB.name(loc);
            final SimpleInterval observed = targetDB.target(name);
            Assert.assertEquals(observed,loc);
            Assert.assertEquals(targetDB.index(name),i);
        }
        Assert.assertNull(targetDB.target("no-interval"));
        Assert.assertEquals(targetDB.index("no-interval"),-1);
    }

    @DataProvider(name="indexRangeData")
    public Object[][] indexRangeData() {
        final List<Object[]> result = new ArrayList<>();

        final Random rdn = new Random(131313);
        for (int i = 0; i < nonOverlappingTargetIntervals.size(); i++) {
            result.add(new Object[]{nonOverlappingTargetIntervals.get(i), new IndexRange(i, i + 1)});
            result.add(new Object[]{
                    TargetsToolsTestUtils.createInterval(nonOverlappingTargetIntervals.get(i).getContig(),
                            nonOverlappingTargetIntervals.get(i).getStart() - 1), new IndexRange(i,i)});
            int j;
            for (j = i + 1; j < nonOverlappingTargetIntervals.size() && rdn.nextBoolean(); j++) {
                if (!nonOverlappingTargetIntervals.get(j).getContig().equals(nonOverlappingTargetIntervals.get(i).getContig())) {
                    break;
                }
            }
            result.add(new Object[]{
                    TargetsToolsTestUtils.createInterval(nonOverlappingTargetIntervals.get(i).getContig(),
                            nonOverlappingTargetIntervals.get(i).getStart(),
                            nonOverlappingTargetIntervals.get(j - 1).getEnd()),
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
                { nonOverlappingTargetIntervals.size() + 1, nonOverlappingTargetIntervals.size() }
        };
    }

    @DataProvider(name="arbitraryTargetRageData")
    public Object[][] arbitraryTargetRageData() {
        final int targetCount = nonOverlappingTargetIntervals.size();
        final Random rdn = new Random(1313); // feeling very lucky.
        final Object[][] result = new Object[100][];
        // make sure there is at least a one element target range request.
        result[0] = new Object[] { targetCount - 4, targetCount - 3 };
        for (int i = 1; i < result.length; i++) {
            final int from = rdn.nextInt(targetCount);
            final int to = rdn.nextInt(targetCount - from) + from;
            result[i] = new Object[] { from , to};
        }
        return result;
    }

    @Test
    public void testRandomLookupTargetByLocation() {
        final Object[][] lookups = targetLookUpData();
        final Random rdn = new Random(11132331);
        for (int i = 0; i < 1000; i++) {
            final Object[] params = lookups[rdn.nextInt(lookups.length)];
            @SuppressWarnings("unchecked")
            final SimpleInterval query = (SimpleInterval) params[0];
            @SuppressWarnings("unchecked")
            final SimpleInterval expected = (SimpleInterval) params[1];
            final SimpleInterval observed = targetDB.target(query);
            Assert.assertEquals(observed,expected);
        }

    }

    @DataProvider(name="targetLookUpData")
    public Object[][] targetLookUpData() {
        final List<Object[]> result = new ArrayList<>();
        for (int i = 0; i < nonOverlappingTargetIntervals.size(); i++) {
            result.add(new Object[]{nonOverlappingTargetIntervals.get(i), nonOverlappingTargetIntervals.get(i), i});
        }
        for (int i = 0; i < nonOverlappingTargetIntervals.size(); i++) {
            final SimpleInterval interval = nonOverlappingTargetIntervals.get(i);
            result.add(new Object[]{new SimpleInterval(interval.getContig(),interval.getStart(),interval.getStart()), interval, i});
        }
        for (int i = 0; i < nonOverlappingTargetIntervals.size(); i++) {
            final SimpleInterval interval = nonOverlappingTargetIntervals.get(i);
            result.add(new Object[]{new SimpleInterval(interval.getContig(),interval.getEnd(),interval.getEnd()), interval, i});
        }
        for (int i = 1; i < nonOverlappingTargetIntervals.size(); i++) {
            final SimpleInterval previous = nonOverlappingTargetIntervals.get(i-1);
            final SimpleInterval next = nonOverlappingTargetIntervals.get(i);
            final SimpleInterval query = previous.getContig().equals(next.getContig())
                    ? TargetsToolsTestUtils.createInterval(previous.getContig(), previous.getEnd() + 1, next.getStart() - 1)
                    : TargetsToolsTestUtils.createInterval(next.getContig(), 1, next.getStart() - 1);
            result.add(new Object[]{query, null, -i - 1});
        }
        return result.toArray(new Object[result.size()][]);
    }

    @DataProvider(name="wrongTargetLookUpData")
    public Object[][] wrongTargetLookUpData() {
        return new Object[][] {
                {null},
                {TargetsToolsTestUtils.createOverEntireContig(nonOverlappingTargetIntervals.get(0).getContig())},
                {TargetsToolsTestUtils.createInterval(nonOverlappingTargetIntervals.get(0).getContig(),
                        nonOverlappingTargetIntervals.get(0).getEnd(),
                        nonOverlappingTargetIntervals.get(1).getStart())}
        };
    }

    @DataProvider(name="trimmedIntervalData")
    public Object[][] trimmedIntervalData() {

        final List<SimpleInterval> intervals = nonOverlappingTargetIntervals;    //alias to reduce clutter

        final String contig = intervals.get(0).getContig();
        //find the last index of the first contig
        int j = 0;
        while (j < intervals.size() && intervals.get(j).getContig().equals(contig)) {
            j++;
        }
        j--;

        //shift of untrimmed interval ends
        final int offset = minimumTargetIntergapSize /2;

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
        Assert.assertEquals(targetDB.trimmedInterval(untrimmed), trimmed);
    }
}
