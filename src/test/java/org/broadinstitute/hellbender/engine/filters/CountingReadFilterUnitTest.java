package org.broadinstitute.hellbender.engine.filters;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public final class CountingReadFilterUnitTest {

    // Mirrors ReadFilterUnitTest
    static final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader(1, 1, 10);
    static final GATKRead goodRead = ArtificialReadUtils.createArtificialRead(header, "Zuul", 0, 2, 2);
    static final GATKRead endBad = ArtificialReadUtils.createArtificialRead(header, "Peter", 0, 1, 100);
    static final GATKRead startBad = ArtificialReadUtils.createArtificialRead(header, "Ray", 0, -1, 2);
    static final GATKRead bothBad = ArtificialReadUtils.createArtificialRead(header, "Egon", 0, -1, 100);

    static final ReadFilter startOk = new ReadFilter() {
        private static final long serialVersionUID = 1L;
        @Override public boolean test(final GATKRead read){return read.getStart() >= 1;}
    };
    static final ReadFilter endOk = new ReadFilter() {
        private static final long serialVersionUID = 1L;
        @Override public boolean test(final GATKRead read){return read.getEnd() <= 10;}
    };

    // Helper to verify post-filtering filter state
    private void verifyFilterState(CountingReadFilter rf, boolean expected) {
        long count = rf.getFilteredCount();
        String rfSummary = rf.getSummaryLine();

        if (expected) {
            Assert.assertTrue(0 == count);
            Assert.assertEquals(0, rfSummary.indexOf("No reads filtered"));
        } else {
            Assert.assertTrue(1 == count);
            Assert.assertEquals(0, rfSummary.indexOf("1 read(s) filtered"));
        }
    }

    @DataProvider(name = "readsStartEnd")
    public Object[][] readsStartEnd() {
        return new Object[][]{
                {goodRead, true, true},
                {startBad, false, true},
                {endBad, true, false},
                {bothBad, false, false}
        };
    }

    @Test(dataProvider = "readsStartEnd")
    public void testTest(GATKRead read, boolean start, boolean end) {

        CountingReadFilter startOkCounting = new CountingReadFilter(startOk);
        Assert.assertEquals(startOkCounting.test(read), start);
        verifyFilterState(startOkCounting, start);

        CountingReadFilter endOkCounting = new CountingReadFilter(endOk);
        Assert.assertEquals(endOkCounting.test(read), end);
        verifyFilterState(endOkCounting, end);
    }

    @Test(dataProvider = "readsStartEnd")
    public void testNegate(GATKRead read, boolean start, boolean end) {
        CountingReadFilter notStartOkCounting = new CountingReadFilter(startOk).negate();
        Assert.assertEquals(notStartOkCounting.test(read), !start);
        verifyFilterState(notStartOkCounting, !start);

        CountingReadFilter notEndOkCounting = new CountingReadFilter(endOk).negate();
        Assert.assertEquals(notEndOkCounting.test(read), !end);
        verifyFilterState(notEndOkCounting, !end);
    }

    @DataProvider(name = "readsAnd")
    public Object[][] readsAnd() {
        return new Object[][]{
                {goodRead, true},
                {startBad, false},
                {endBad, false},
                {bothBad, false}
        };
    }

    @Test(dataProvider = "readsAnd")
    public void testAnd(GATKRead read, boolean expected) {

        CountingReadFilter startAndEndOk = new CountingReadFilter(startOk).and(new CountingReadFilter(endOk));
        Assert.assertEquals(startAndEndOk.test(read), expected);
        verifyFilterState(startAndEndOk, expected);

        CountingReadFilter endAndStartOk = new CountingReadFilter(endOk).and(new CountingReadFilter(startOk));
        Assert.assertEquals(endAndStartOk.test(read), expected);
        verifyFilterState(endAndStartOk, expected);
    }

    @DataProvider(name = "readsOr")
    public Object[][] readsOr() {
        return new Object[][]{
                {goodRead, true},
                {startBad, true},
                {endBad, true},
                {bothBad, false}
        };
    }

    @Test(dataProvider = "readsOr")
    public void testOr(GATKRead read, boolean expected) {

        CountingReadFilter startOrEndOk = new CountingReadFilter(startOk).or(new CountingReadFilter(endOk));
        Assert.assertEquals(startOrEndOk.test(read), expected);
        verifyFilterState(startOrEndOk, expected);

        CountingReadFilter endOrStartOk = new CountingReadFilter(endOk).or(new CountingReadFilter(startOk));
        Assert.assertEquals(endOrStartOk.test(read), expected);
        verifyFilterState(endOrStartOk, expected);
    }

    @DataProvider(name = "deeper")
    public Object[][] deeper() {
        return new Object[][]{
                {goodRead, false},
                {startBad, true},
                {endBad, true},
                {bothBad, false}
        };
    }

    private CountingReadFilter readChecksOut() {
        return new CountingReadFilter(startOk)
                .or(new CountingReadFilter(endOk))
                .and(new CountingReadFilter(new ReadFilter() {
                    private static final long serialVersionUID = 1L;
                    @Override public boolean test(final GATKRead read){return !read.getName().equals("Zuul");}
                }));
    }

    @Test(dataProvider = "deeper")
    public void testDeeperChaining(GATKRead read, boolean expected) {

        CountingReadFilter readCheckOutCounting = readChecksOut();
        Assert.assertEquals(readCheckOutCounting.test(read), expected);
        verifyFilterState(readCheckOutCounting, expected);

        readCheckOutCounting = readChecksOut().and(readChecksOut());
        Assert.assertEquals(readCheckOutCounting.test(read), expected);
        verifyFilterState(readCheckOutCounting, expected);

        readCheckOutCounting = readChecksOut().and(new CountingReadFilter(
                new ReadFilter() {
                    private static final long serialVersionUID = 1L;
                    @Override public boolean test(final GATKRead read){return false;}
                }));
        Assert.assertEquals(readCheckOutCounting.test(read), false);
        verifyFilterState(readCheckOutCounting, false);

        readCheckOutCounting = readChecksOut().or(new CountingReadFilter(
                new ReadFilter() {
                    private static final long serialVersionUID = 1L;
                    @Override public boolean test(final GATKRead read){return true;}
                }));
        Assert.assertEquals(readCheckOutCounting.test(read), true);
        verifyFilterState(readCheckOutCounting, true);
    }

    @DataProvider(name = "multipleRejection")
    public Object[][] multipleRejection() {
        return new Object[][] {
            {
                new GATKRead[] { goodRead, goodRead, goodRead }, 0
            },
            {
                new GATKRead[] { goodRead, goodRead, bothBad }, 1
            },
            {
                new GATKRead[] { startBad, bothBad, goodRead, startBad }, 3
            },
        };
    }

    @Test(dataProvider = "multipleRejection")
    public void testRootFilterCounts(GATKRead[] reads, int expectedRejectionCount) {

        CountingReadFilter startOkCounting = new CountingReadFilter(startOk);
        Arrays.asList(reads).stream().filter(startOkCounting).count();  // force the stream to be consumed
        Assert.assertEquals(startOkCounting.getFilteredCount(), expectedRejectionCount);
        startOkCounting.resetFilteredCount();
        Assert.assertEquals(startOkCounting.getFilteredCount(), 0);
    }

    @DataProvider(name = "subFilterCounts")
    public Object[][] subFilterCounts() {
        return new Object[][] {
                {
                        new GATKRead[]{goodRead, startBad, bothBad, bothBad}, 1L, 2L, 1L
                },
                {
                        new GATKRead[]{goodRead, goodRead, goodRead, bothBad }, 3L, 3L, 3L
                },
                {
                        new GATKRead[]{goodRead, startBad, endBad, bothBad}, 2L, 3L, 2L
                },
        };
    }

    @Test(dataProvider = "subFilterCounts")
    public void testSubFilterCounts(GATKRead[] reads, long totalRejections, long startEndRejections, long nameRejections) {

        CountingReadFilter badStart = new CountingReadFilter(
                new ReadFilter() {
                    private static final long serialVersionUID = 1L;
                    @Override public boolean test(final GATKRead read){return read.getStart() < 1;}
                }
        );
        CountingReadFilter badEnd = new CountingReadFilter(
                new ReadFilter() {
                    private static final long serialVersionUID = 1L;
                    @Override public boolean test(final GATKRead read){return read.getEnd() > 10;}
                }
        );
        CountingReadFilter badStartAndEnd = badStart.and(badEnd);

        CountingReadFilter isRay= new CountingReadFilter(
                new ReadFilter() {
                    private static final long serialVersionUID = 1L;
                    @Override public boolean test(final GATKRead read){return read.getName().equals("Ray");}
                }
        );
        CountingReadFilter isEgon = new CountingReadFilter(
                new ReadFilter() {
                    private static final long serialVersionUID = 1L;
                    @Override public boolean test(final GATKRead read){return read.getName().equals("Egon");}
                }
        );
        CountingReadFilter isRayOrEgon = isRay.or(isEgon);

        CountingReadFilter compoundFilter = badStartAndEnd.or(isRayOrEgon);

        Arrays.asList(reads).stream().filter(compoundFilter).count(); // force the stream to be consumed

        Assert.assertEquals(compoundFilter.getFilteredCount(), totalRejections);
        Assert.assertEquals(badStartAndEnd.getFilteredCount(), startEndRejections);
        Assert.assertEquals(isRayOrEgon.getFilteredCount(), nameRejections);

        // test if reset filtered count is correctly propagated
        compoundFilter.resetFilteredCount();
        Assert.assertEquals(compoundFilter.getFilteredCount(), 0);
        Assert.assertEquals(badStartAndEnd.getFilteredCount(), 0);
        Assert.assertEquals(isRayOrEgon.getFilteredCount(), 0);
        Assert.assertEquals(badStart.getFilteredCount(), 0);
        Assert.assertEquals(badEnd.getFilteredCount(), 0);
        Assert.assertEquals(isRay.getFilteredCount(), 0);
        Assert.assertEquals(isEgon.getFilteredCount(), 0);
    }

    @Test
    public void testFromListNull() {
        CountingReadFilter rf = CountingReadFilter.fromList(null, ArtificialReadUtils.createArtificialSamHeader(1, 1, 10));
        Assert.assertTrue(rf.delegateFilter.getClass() == ReadFilterLibrary.AllowAllReadsReadFilter.class);
    }
    @Test
    public void testFromListEmpty() {
        CountingReadFilter rf = CountingReadFilter.fromList(Collections.emptyList(), ArtificialReadUtils.createArtificialSamHeader(1, 1, 10));
        Assert.assertTrue(rf.delegateFilter.getClass() == ReadFilterLibrary.AllowAllReadsReadFilter.class);
    }

    @Test
    public void testFromListSingle() {
        List<ReadFilter> filters = new ArrayList<>();
        filters.add(ReadFilterLibrary.MAPPED);
        CountingReadFilter rf = CountingReadFilter.fromList(filters, ArtificialReadUtils.createArtificialSamHeader(1, 1, 10));
        Assert.assertTrue(rf.delegateFilter.getClass() == ReadFilterLibrary.MAPPED.getClass());
    }

    @Test
    public void testFromListMultiOrdered() {
        List<ReadFilter> filters = new ArrayList<>();
        filters.add(ReadFilterLibrary.MAPPING_QUALITY_AVAILABLE);
        filters.add(ReadFilterLibrary.MAPPED);
        filters.add(ReadFilterLibrary.GOOD_CIGAR);

        // Since we want to ensure that order of the input is honored, we need to test the
        // structure of the filter rather than the result
        CountingReadFilter rf = CountingReadFilter.fromList(filters, ArtificialReadUtils.createArtificialSamHeader(1, 1, 10));

        Assert.assertTrue(rf.getClass() == CountingReadFilter.CountingAndReadFilter.class);
        CountingReadFilter.CountingAndReadFilter andFilter = (CountingReadFilter.CountingAndReadFilter) rf;

        // lhs is a Counting and filter; rhs is a counting filter that delegates to GOOD_CIGAR
        Assert.assertTrue(andFilter.lhs.getClass() == CountingReadFilter.CountingAndReadFilter.class);
        Assert.assertTrue(andFilter.rhs.delegateFilter.getClass() == ReadFilterLibrary.GOOD_CIGAR.getClass());
        andFilter = (CountingReadFilter.CountingAndReadFilter) andFilter.lhs;

        // lhs is a Counting filter that delegates to MAPPING_QUALITY_AVAILABLE; rhs is a
        // counting filter that delegates to MAPPED
        Assert.assertTrue(andFilter.lhs.getClass() == CountingReadFilter.class);
        Assert.assertTrue(andFilter.lhs.delegateFilter.getClass() == ReadFilterLibrary.MAPPING_QUALITY_AVAILABLE.getClass());
        Assert.assertTrue(andFilter.rhs.getClass() == CountingReadFilter.class);
        Assert.assertTrue(andFilter.rhs.delegateFilter.getClass() == ReadFilterLibrary.MAPPED.getClass());
    }

}

