package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.utils.IndexRange;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import org.testng.internal.Utils;

import java.lang.reflect.Method;
import java.util.*;

/**
 * Unit tests for {@link ExonCollection}'s default method implementations.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class ExonCollectionUnitTest extends BaseTest {

    /**
     * General exon_db_test_data, will contain ExonDB instances as a single arguments.
     *
     *  <p>
     *  Initialized by {@link #setUp} before any tests.
     * </p>
     */
    private static Object[][] EXON_DB_TEST_DATA;

    @BeforeClass
    public void setUp() {
        final Random rdn = new Random(11131719);
        final List<Object[]> exonDBTestData = new ArrayList<>(30);
        exonDBTestData.add(new Object[] { new ExonCollectionStub(0,rdn) });
        exonDBTestData.add(new Object[] { new ExonCollectionStub(1,rdn) });
        exonDBTestData.add(new Object[] { new ExonCollectionStub(33,rdn) });
        final int maxExonCount = ExomeToolsTestUtils.REFERENCE_DICTIONARY.getSequence(0).getSequenceLength() / 50;
        final int minExonCount = 2;
        while (exonDBTestData.size() < 30) {
            final int exonCount = rdn.nextInt(maxExonCount - minExonCount + 1) + minExonCount;
            exonDBTestData.add(new Object[] { new ExonCollectionStub(exonCount,rdn)});
        }
        EXON_DB_TEST_DATA = exonDBTestData.toArray(new Object[exonDBTestData.size()][]);
    }

    @Test(dataProvider = "exonDBData")
    public void testExonByName(final ExonCollection<SimpleInterval> exonCollection) {
        Assert.assertNull(exonCollection.exon("silly-name"));
        for (int i = 0; i < exonCollection.exonCount(); i++) {
            Assert.assertEquals(exonCollection.exon("" + i), exonCollection.exon(i));
        }
    }

    @Test(dataProvider = "exonDBData")
    public void testExonByLocationExactMatch(final ExonCollection<SimpleInterval> exonCollection) {
        for (int i = 0; i < exonCollection.exonCount(); i++) {
            final SimpleInterval si = exonCollection.exon(i);
            Assert.assertEquals(exonCollection.exon(si), si);
        }
    }

    @Test(dataProvider = "exonDBData")
    public void testExonByLocationNonOverlapping(final ExonCollection<SimpleInterval> exonCollection) {
        final SimpleInterval otherChromosome =  ExomeToolsTestUtils.createOverEntireContig(2);
        Assert.assertNull(exonCollection.exon(otherChromosome));
        for (int i = 0; i < exonCollection.exonCount() - 1; i++) {
            final int start = exonCollection.exon(i).getEnd() + 1;
            final int end = exonCollection.exon(i+1).getStart() - 1;
            if (end < start) {   // avoid zero-length intra-exon intervals (due to artificial data construction).
                continue;
            }
            final SimpleInterval loc = ExomeToolsTestUtils.createInterval(exonCollection.exon(i).getContig(), start, end);
            Assert.assertNull(exonCollection.exon(loc),"query = " + exonCollection.exon(i).getContig() + ":" + start + "-" + end);
        }
    }

    @Test(dataProvider = "twoOrMoreExonDBData")
    public void testExonByLocationAmbiguousOverlapping(final ExonCollection<SimpleInterval> exonCollection) {
        for (int i = 0; i < exonCollection.exonCount() - 1; i++) {
            final int start = exonCollection.exon(i).getEnd() - 1;
            final int end = exonCollection.exon(i+1).getStart() + 1;
            final SimpleInterval loc = ExomeToolsTestUtils.createInterval(exonCollection.exon(i).getContig(), start, end);
            try {
                exonCollection.exon(loc);
                Assert.fail("expected exception. Index == " + i);
            } catch (final ExonCollection.AmbiguousExonException ex) {
                // ok.
            } catch (final Throwable t) {
                Assert.fail("wrong exception: ",t);
            }
        }
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testForEachNullTask() {
        final ExonCollection<SimpleInterval> exonCollection = new ExonCollectionStub(0,new Random(13));
        exonCollection.forEachExon(ExomeToolsTestUtils.createOverEntireContig(0), null);
    }

    @Test(dataProvider = "exonDBData")
    public void testForEachSingleCallExactMatch(final ExonCollection<SimpleInterval> exonCollection) {
        for (int i = 0; i < exonCollection.exonCount(); i++) {
            final SimpleInterval simpleInterval = exonCollection.exon(i);
            final SimpleInterval loc = ExomeToolsTestUtils.createInterval(simpleInterval.getContig(), simpleInterval.getStart(), simpleInterval.getEnd());
            final int expectedIdx = i;
            final int[] totalCalls = new int[] { 0 };
            exonCollection.forEachExon(loc,(idx,exon) -> {
                Assert.assertEquals(exon,simpleInterval);
                Assert.assertEquals(idx,expectedIdx);
                totalCalls[0]++;
            });
            Assert.assertEquals(totalCalls[0],1);
        }
    }

    @Test(dataProvider = "exonDBData")
    public void testForEachExhaustiveRanges(final ExonCollection<SimpleInterval> exonCollection) {
        final int min = 0;
        final int max = exonCollection.exonCount();
        for (int i = min; i < max; i++) {

            final SimpleInterval iInterval = exonCollection.exon(i);
            final String contig = iInterval.getContig();
            final int start = (iInterval.getStart() + iInterval.getEnd()) >> 1;
            for (int j = i; j < max; j++) {
                final SimpleInterval jInterval = exonCollection.exon(j);
                if (!jInterval.getContig().equals(contig)) {
                    break;
                }
                final int end = (jInterval.getStart() + jInterval.getEnd() + 1) >> 1;
                if (start > end) {
                    continue;
                }
                final SimpleInterval loc = ExomeToolsTestUtils.createInterval(contig, start, end);
                final List<SimpleInterval> visited = new ArrayList<>(j - i + 1);
                final int iFinal = i;
                exonCollection.forEachExon(loc, (idx, e) -> {
                    Assert.assertEquals(idx, iFinal + visited.size());
                    Assert.assertEquals(e, exonCollection.exon(idx));
                    visited.add(e);
                });
                Assert.assertEquals(visited.size(),j - i + 1);
            }
        }
    }

    @Test(dataProvider = "exonDBData", expectedExceptions = IllegalArgumentException.class)
    public void testExonsByBadIndexRange(final ExonCollection<SimpleInterval> exonCollection) {
        exonCollection.exons(new IndexRange(exonCollection.exonCount(), exonCollection.exonCount() + 1));
    }

    @Test(dataProvider = "exonDBData", expectedExceptions = IllegalArgumentException.class)
    public void testExonsByNullIndexRange(final ExonCollection<SimpleInterval> exonCollection) {
        exonCollection.exons((IndexRange)null);
    }


    @Test(dataProvider = "exonDBData")
    public void testForEachExhaustiveZeroRanges(final ExonCollection<SimpleInterval> exonCollection) {
        final int min = 0;
        final int max = exonCollection.exonCount();
        for (int i = min; i < max - 1; i++) {
            final SimpleInterval iInterval = exonCollection.exon(i);
            final SimpleInterval jInterval = exonCollection.exon(i + 1);
            if (!jInterval.getContig().equals(iInterval.getContig()))
                continue;
            final int start = iInterval.getEnd() + 1;
            final int end = jInterval.getStart() - 1;
            if (start < end)
                continue;
            exonCollection.forEachExon(
                    ExomeToolsTestUtils.createInterval(iInterval.getContig(), start, end),
                    (idx, e) -> Assert.fail("no exon should be call"));
        }
    }

    @Test(dataProvider = "exonDBData")
    public void testExonsByLocationExhaustiveZeroRanges(final ExonCollection<SimpleInterval> exonCollection) {
        final int min = 0;
        final int max = exonCollection.exonCount();
        for (int i = min; i < max - 1; i++) {
            final SimpleInterval iInterval = exonCollection.exon(i);
            final SimpleInterval jInterval = exonCollection.exon(i + 1);
            if (!jInterval.getContig().equals(iInterval.getContig()))
                continue;
            final int start = iInterval.getEnd() + 1;
            final int end = jInterval.getStart() - 1;
            if (start < end)
                continue;
            final List<SimpleInterval> observed = exonCollection.exons(ExomeToolsTestUtils.createInterval(iInterval.getContig(), start, end));
            Assert.assertEquals(observed.size(), 0);
        }
    }

    @Test(dataProvider = "exonDBData")
    public void testExonSize(final ExonCollection<SimpleInterval> exonCollection) {
        long totalSize = 0;
        for (int i = 0; i < exonCollection.exonCount(); i++) {
            totalSize += exonCollection.exon(i).size();
        }
        Assert.assertEquals(exonCollection.exomeSize(), totalSize);
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
     * Stub class to test default implementations of {@link ExonCollection}.
     *
     * <p>
     *     This class should not implement any method that has a default implementation in {@code ExonDB} as
     *     that would interfere with its testing.
     * </p>
     */
    public static class ExonCollectionStub implements ExonCollection<SimpleInterval> {

        private final List<SimpleInterval> intervals;

        public ExonCollectionStub(final int numberOfExons, final Random rdn) {
            if (numberOfExons == 0) {
                this.intervals = Collections.emptyList();
                return;
            }
            List<SimpleInterval> intervals = new ArrayList<>(numberOfExons);
            final String contig = ExomeToolsTestUtils.REFERENCE_DICTIONARY.getSequence(0).getSequenceName();
            final int contigLength = ExomeToolsTestUtils.REFERENCE_DICTIONARY.getSequence(0).getSequenceLength();
            final float avgSlotSize = contigLength / numberOfExons;

            int nextBasePos = 0;
            for (int i = 0; i < numberOfExons; i++) {
                final int slotLength = (int) Math.max(1,Math.min(3 * avgSlotSize,(rdn.nextGaussian() * avgSlotSize * .5) + avgSlotSize));
                final int exonLength = (int) Math.max(1,slotLength * rdn.nextDouble());
                final int start = Math.max(nextBasePos,1 + nextBasePos + ((slotLength - exonLength) >> 1));
                final int stop = start + exonLength;
                intervals.add(new SimpleInterval(contig,start,stop));
                nextBasePos += slotLength;
                nextBasePos = Math.max(nextBasePos,stop + 1);
            }
            this.intervals = Collections.unmodifiableList(intervals);
        }

        @Override
        public int exonCount() {
            return intervals.size();
        }

        @Override
        public SimpleInterval exon(final int index) {
            return intervals.get(index);
        }

        @Override
        public String name(final SimpleInterval exon) {
            if (exon == null) {
                throw new IllegalArgumentException("the input exon cannot be null");
            }
            return exon.toString();
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
        public List<SimpleInterval> exons() {
            return intervals;
        }

        @Override
        public SimpleInterval location(final int index) {
            if (index < 0) {
                throw new IndexOutOfBoundsException();
            } if (index >= exonCount()) {
                throw new IndexOutOfBoundsException();
            }
            return intervals.get(index);
        }

        @Override
        public SimpleInterval location(final SimpleInterval exon) {
            return exon;
        }

        @Override
        public IndexRange indexRange(final SimpleInterval location) {
            if (!location.getContig().equals(ExomeToolsTestUtils.REFERENCE_DICTIONARY.getSequence(0).getSequenceName()))
                return new IndexRange(exonCount(),exonCount());
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
    public void testExonDBSubOnlyImplementsNonDefaults() throws NoSuchMethodException {
        final List<Method> overriddenDefaults = new ArrayList<>(10);
        for (final Method method : ExonCollection.class.getMethods()) {
            if (!method.isDefault())
                continue;
            final Method implementationMethod = ExonCollectionStub.class.getMethod(method.getName(), method.getParameterTypes());
            if (!implementationMethod.getDeclaringClass().isAssignableFrom(ExonCollection.class))
                overriddenDefaults.add(method);
        }
        Assert.assertTrue(overriddenDefaults.isEmpty(), "Bug in testing code!!! " + ExonCollectionStub.class.getName() +
                " must not override any default method in " + ExonCollection.class + " as we present to test those. " +
                "Currently it overrides the following default methods: "
                + Utils.join(Arrays.asList(overriddenDefaults.stream().map(Method::getName).toArray()), ", "));
    }

    @DataProvider(name="exonDBData")
    public Object[][] exonDBData() {
        return EXON_DB_TEST_DATA;
    }

    @DataProvider(name="twoOrMoreExonDBData")
    public Object[][] twoOrMoreExonDBData() {
        final List<Object[]> result = new ArrayList<>(EXON_DB_TEST_DATA.length);
        for (final Object[] params : EXON_DB_TEST_DATA) {
            @SuppressWarnings("unchecked")
            final ExonCollection<SimpleInterval> db = (ExonCollection<SimpleInterval>)params[0];
            if (db.exonCount() <= 1)
                continue;
            result.add(params);
        }
        return result.toArray(new Object[result.size()][]);
    }
}
