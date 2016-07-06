package org.broadinstitute.hellbender.tools.exome;

import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.utils.IndexRange;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.test.TargetsToolsTestUtils;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.lang.reflect.Method;
import java.util.*;

/**
 * Unit tests for {@link TargetCollection}'s default method implementations.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class TargetCollectionUnitTest extends BaseTest {

    /**
     * General target_db_test_data, will contain TargetDB instances as a single arguments.
     *
     *  <p>
     *  Initialized by {@link #setUp} before any tests.
     * </p>
     */
    private static Object[][] TARGET_DB_TEST_DATA;

    @BeforeClass
    public void setUp() {
        final Random rdn = new Random(11131719);
        final List<Object[]> targetDBTestData = new ArrayList<>(30);
        targetDBTestData.add(new Object[]{new TargetCollectionStub(0, rdn)});
        targetDBTestData.add(new Object[] { new TargetCollectionStub(1,rdn) });
        targetDBTestData.add(new Object[] { new TargetCollectionStub(33,rdn) });
        final int maxTargetCount = TargetsToolsTestUtils.REFERENCE_DICTIONARY.getSequence(0).getSequenceLength() / 50;
        final int minTargetCount = 2;
        while (targetDBTestData.size() < 30) {
            final int targetCount = rdn.nextInt(maxTargetCount - minTargetCount + 1) + minTargetCount;
            targetDBTestData.add(new Object[] { new TargetCollectionStub(targetCount,rdn)});
        }
        TARGET_DB_TEST_DATA = targetDBTestData.toArray(new Object[targetDBTestData.size()][]);
    }

    @Test(dataProvider = "targetDBData")
    public void testTargetByName(final TargetCollection<SimpleInterval> targetCollection) {
        Assert.assertNull(targetCollection.target("silly-name"));
        for (int i = 0; i < targetCollection.targetCount(); i++) {
            Assert.assertEquals(targetCollection.target("" + i), targetCollection.target(i));
        }
    }

    @Test(dataProvider = "targetDBData")
    public void testTargetByLocationExactMatch(final TargetCollection<SimpleInterval> targetCollection) {
        for (int i = 0; i < targetCollection.targetCount(); i++) {
            final SimpleInterval si = targetCollection.target(i);
            Assert.assertEquals(targetCollection.target(si), si);
        }
    }

    @Test(dataProvider = "targetDBData")
    public void testTargetByLocationNonOverlapping(final TargetCollection<SimpleInterval> targetCollection) {
        final SimpleInterval otherChromosome =  TargetsToolsTestUtils.createOverEntireContig(2);
        Assert.assertNull(targetCollection.target(otherChromosome));
        for (int i = 0; i < targetCollection.targetCount() - 1; i++) {
            final int start = targetCollection.target(i).getEnd() + 1;
            final int end = targetCollection.target(i + 1).getStart() - 1;
            if (end < start) {   // avoid zero-length intra-target intervals (due to artificial data construction).
                continue;
            }
            final SimpleInterval loc = TargetsToolsTestUtils.createInterval(targetCollection.target(i).getContig(), start, end);
            Assert.assertNull(targetCollection.target(loc),"query = " + targetCollection.target(i).getContig() + ":" + start + "-" + end);
        }
    }

    @Test(dataProvider = "twoOrMoreTargetDBData")
    public void testTargetByLocationAmbiguousOverlapping(final TargetCollection<SimpleInterval> targetCollection) {
        for (int i = 0; i < targetCollection.targetCount() - 1; i++) {
            final int start = targetCollection.target(i).getEnd() - 1;
            final int end = targetCollection.target(i + 1).getStart() + 1;
            final SimpleInterval loc = TargetsToolsTestUtils.createInterval(targetCollection.target(i).getContig(), start, end);
            try {
                targetCollection.target(loc);
                Assert.fail("expected exception. Index == " + i);
            } catch (final TargetCollection.AmbiguousTargetException ex) {
                // ok.
            } catch (final Throwable t) {
                Assert.fail("wrong exception: ",t);
            }
        }
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testForEachNullTask() {
        final TargetCollection<SimpleInterval> targetCollection = new TargetCollectionStub(0,new Random(13));
        targetCollection.forEachTarget(TargetsToolsTestUtils.createOverEntireContig(0), null);
    }

    @Test(dataProvider = "targetDBData")
    public void testForEachSingleCallExactMatch(final TargetCollection<SimpleInterval> targetCollection) {
        for (int i = 0; i < targetCollection.targetCount(); i++) {
            final SimpleInterval simpleInterval = targetCollection.target(i);
            final SimpleInterval loc = TargetsToolsTestUtils.createInterval(simpleInterval.getContig(), simpleInterval.getStart(), simpleInterval.getEnd());
            final int expectedIdx = i;
            final int[] totalCalls = new int[] { 0 };
            targetCollection.forEachTarget(loc, (idx, target) -> {
                Assert.assertEquals(target, simpleInterval);
                Assert.assertEquals(idx, expectedIdx);
                totalCalls[0]++;
            });
            Assert.assertEquals(totalCalls[0],1);
        }
    }

    @Test(dataProvider = "targetDBData")
    public void testForEachExhaustiveRanges(final TargetCollection<SimpleInterval> targetCollection) {
        final int min = 0;
        final int max = targetCollection.targetCount();
        for (int i = min; i < max; i++) {

            final SimpleInterval iInterval = targetCollection.target(i);
            final String contig = iInterval.getContig();
            final int start = (iInterval.getStart() + iInterval.getEnd()) >> 1;
            for (int j = i; j < max; j++) {
                final SimpleInterval jInterval = targetCollection.target(j);
                if (!jInterval.getContig().equals(contig)) {
                    break;
                }
                final int end = (jInterval.getStart() + jInterval.getEnd() + 1) >> 1;
                if (start > end) {
                    continue;
                }
                final SimpleInterval loc = TargetsToolsTestUtils.createInterval(contig, start, end);
                final List<SimpleInterval> visited = new ArrayList<>(j - i + 1);
                final int iFinal = i;
                targetCollection.forEachTarget(loc, (idx, e) -> {
                    Assert.assertEquals(idx, iFinal + visited.size());
                    Assert.assertEquals(e, targetCollection.target(idx));
                    visited.add(e);
                });
                Assert.assertEquals(visited.size(),j - i + 1);
            }
        }
    }

    @Test(dataProvider = "targetDBData", expectedExceptions = IllegalArgumentException.class)
    public void testTargetsByBadIndexRange(final TargetCollection<SimpleInterval> targetCollection) {
        targetCollection.targets(new IndexRange(targetCollection.targetCount(), targetCollection.targetCount() + 1));
    }

    @Test(dataProvider = "targetDBData", expectedExceptions = IllegalArgumentException.class)
    public void testTargetsByNullIndexRange(final TargetCollection<SimpleInterval> targetCollection) {
        targetCollection.targets((IndexRange) null);
    }


    @Test(dataProvider = "targetDBData")
    public void testForEachExhaustiveZeroRanges(final TargetCollection<SimpleInterval> targetCollection) {
        final int min = 0;
        final int max = targetCollection.targetCount();
        for (int i = min; i < max - 1; i++) {
            final SimpleInterval iInterval = targetCollection.target(i);
            final SimpleInterval jInterval = targetCollection.target(i + 1);
            if (!jInterval.getContig().equals(iInterval.getContig()))
                continue;
            final int start = iInterval.getEnd() + 1;
            final int end = jInterval.getStart() - 1;
            if (start < end)
                continue;
            targetCollection.forEachTarget(
                    TargetsToolsTestUtils.createInterval(iInterval.getContig(), start, end),
                    (idx, e) -> Assert.fail("no target should be call"));
        }
    }

    @Test(dataProvider = "targetDBData")
    public void testTargetsByLocationExhaustiveZeroRanges(final TargetCollection<SimpleInterval> targetCollection) {
        final int min = 0;
        final int max = targetCollection.targetCount();
        for (int i = min; i < max - 1; i++) {
            final SimpleInterval iInterval = targetCollection.target(i);
            final SimpleInterval jInterval = targetCollection.target(i + 1);
            if (!jInterval.getContig().equals(iInterval.getContig()))
                continue;
            final int start = iInterval.getEnd() + 1;
            final int end = jInterval.getStart() - 1;
            if (start < end)
                continue;
            final List<SimpleInterval> observed = targetCollection.targets(TargetsToolsTestUtils.createInterval(iInterval.getContig(), start, end));
            Assert.assertEquals(observed.size(), 0);
        }
    }

    @Test(dataProvider = "targetDBData")
    public void testTargetSize(final TargetCollection<SimpleInterval> targetCollection) {
        long totalSize = 0;
        for (int i = 0; i < targetCollection.targetCount(); i++) {
            totalSize += targetCollection.target(i).size();
        }
        Assert.assertEquals(targetCollection.totalSize(), totalSize);
    }

    @Test
    public void testIndexRangeCreation() {
       for (int i = 0; i < 10; i++)
           for (int j = i; j < 10; j++) {
               final IndexRange range = new IndexRange(i, j);
               Assert.assertEquals(range.from, i);
               Assert.assertEquals(range.to, j);
           }
    }

    @Test
    public void testIndexRangeSize() {
        for (int i = 0; i < 10; i++)
            for (int j = i; j < 10; j++) {
                final IndexRange range = new IndexRange(i, j);
                Assert.assertEquals(range.size(),j - i);
            }
    }

    @Test
    public void testIndexRangeEqual() {
        final IndexRange r1 = new IndexRange(1,2);
        final IndexRange r1_bis = new IndexRange(1,2);
        final IndexRange r2 = new IndexRange(2,2);
        final IndexRange r3 = new IndexRange(1,10);
        Assert.assertEquals(r1,r1_bis);
        Assert.assertEquals(r1,r1);
        Assert.assertNotEquals(r1,r2);
        Assert.assertNotEquals(r1,null);
        Assert.assertNotEquals(r1, new Object());
        Assert.assertNotEquals(r1, r3);
    }

    @Test
    public void testIndexRangeHashCode() {
        final IndexRange r1 = new IndexRange(1,2);
        final IndexRange r1_bis = new IndexRange(1,2);
        Assert.assertEquals(r1.hashCode(),r1_bis.hashCode());
        Assert.assertEquals(r1.hashCode(),r1.hashCode());
    }

    /**
     * Stub class to test default implementations of {@link TargetCollection}.
     *
     * <p>
     *     This class should not implement any method that has a default implementation in {@link TargetCollection} as
     *     that would interfere with its testing.
     * </p>
     */
    public static class TargetCollectionStub implements TargetCollection<SimpleInterval> {

        private final List<SimpleInterval> intervals;

        public TargetCollectionStub(final int numberOfTargets, final Random rdn) {
            if (numberOfTargets == 0) {
                this.intervals = Collections.emptyList();
                return;
            }
            List<SimpleInterval> intervals = new ArrayList<>(numberOfTargets);
            final String contig = TargetsToolsTestUtils.REFERENCE_DICTIONARY.getSequence(0).getSequenceName();
            final int contigLength = TargetsToolsTestUtils.REFERENCE_DICTIONARY.getSequence(0).getSequenceLength();
            final float avgSlotSize = (float) contigLength / numberOfTargets;

            int nextBasePos = 0;
            for (int i = 0; i < numberOfTargets; i++) {
                final int slotLength = (int) Math.max(1,Math.min(3 * avgSlotSize,(rdn.nextGaussian() * avgSlotSize * .5) + avgSlotSize));
                final int targetLength = (int) Math.max(1,slotLength * rdn.nextDouble());
                final int start = Math.max(nextBasePos,1 + nextBasePos + ((slotLength - targetLength) >> 1));
                final int stop = start + targetLength;
                intervals.add(new SimpleInterval(contig,start,stop));
                nextBasePos += slotLength;
                nextBasePos = Math.max(nextBasePos,stop + 1);
            }
            this.intervals = Collections.unmodifiableList(intervals);
        }

        @Override
        public int targetCount() {
            return intervals.size();
        }

        @Override
        public SimpleInterval target(final int index) {
            return intervals.get(index);
        }

        @Override
        public String name(final SimpleInterval target) {
            Utils.nonNull(target, "the input target cannot be null");
            return target.toString();
        }

        @Override
        public int index(final String name) {
            if (name.matches("^\\d+$")) {
                int index = Integer.valueOf(name);
                if (index >= intervals.size()) {
                    return -1;
                } else {
                    return index;
                }
            } else {
                return -1;
            }
        }

        @Override
        public List<SimpleInterval> targets() {
            return intervals;
        }

        @Override
        public SimpleInterval location(final int index) {
            if (index < 0) {
                throw new IndexOutOfBoundsException();
            } if (index >= targetCount()) {
                throw new IndexOutOfBoundsException();
            }
            return intervals.get(index);
        }

        @Override
        public SimpleInterval location(final SimpleInterval target) {
            return Utils.nonNull(target, "the input target cannot be null");
        }

        @Override
        public IndexRange indexRange(final Locatable location) {
            if (!location.getContig().equals(TargetsToolsTestUtils.REFERENCE_DICTIONARY.getSequence(0).getSequenceName()))
                return new IndexRange(targetCount(), targetCount());
            else {
                int from = 0;
                while (from < intervals.size()) {
                    final SimpleInterval nextInterval = intervals.get(from);
                    if (nextInterval.getEnd() >= location.getStart()) {
                        break;
                    }
                    from++;
                }
                int to = from;
                while (to < intervals.size()) {
                    final SimpleInterval nextInterval = intervals.get(to);
                    if (nextInterval.getStart() > location.getEnd()) {
                        break;
                    }
                    to++;
                }
                return new IndexRange(from,to);
            }
        }
    }

    @Test
    public void testTargetDBSubOnlyImplementsNonDefaults() throws NoSuchMethodException {
        final List<Method> overriddenDefaults = new ArrayList<>(10);
        for (final Method method : TargetCollection.class.getMethods()) {
            if (!method.isDefault())
                continue;
            final Method implementationMethod = TargetCollectionStub.class.getMethod(method.getName(), method.getParameterTypes());
            if (!implementationMethod.getDeclaringClass().isAssignableFrom(TargetCollection.class))
                overriddenDefaults.add(method);
        }
        Assert.assertTrue(overriddenDefaults.isEmpty(), "Bug in testing code!!! " + TargetCollectionStub.class.getName() +
                " must not override any default method in " + TargetCollection.class + " as we present to test those. " +
                "Currently it overrides the following default methods: "
                + Utils.join(", ",Arrays.asList(overriddenDefaults.stream().map(Method::getName).toArray())));
    }

    @DataProvider(name="targetDBData")
    public Object[][] targetDBData() {
        return TARGET_DB_TEST_DATA;
    }

    @DataProvider(name="twoOrMoreTargetDBData")
    public Object[][] twoOrMoreTargetDBData() {
        final List<Object[]> result = new ArrayList<>(TARGET_DB_TEST_DATA.length);
        for (final Object[] params : TARGET_DB_TEST_DATA) {
            @SuppressWarnings("unchecked")
            final TargetCollection<SimpleInterval> db = (TargetCollection<SimpleInterval>)params[0];
            if (db.targetCount() <= 1)
                continue;
            result.add(params);
        }
        return result.toArray(new Object[result.size()][]);
    }
}
