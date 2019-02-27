package org.broadinstitute.hellbender.engine.filters;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public final class ReadFilterUnitTest {

    // Mirrors CountingReadFilterUnitTest
    static final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader(1, 1, 10);
    static final GATKRead goodRead = ArtificialReadUtils.createArtificialRead(header, "Zuul", 0, 2,2);
    static final GATKRead endBad = ArtificialReadUtils.createArtificialRead(header, "Peter", 0, 1,100);
    static final GATKRead startBad = ArtificialReadUtils.createArtificialRead(header, "Ray", 0, -1,2);
    static final GATKRead bothBad = ArtificialReadUtils.createArtificialRead(header, "Egon", 0, -1,100);
    static final ReadFilter startOk = new ReadFilter() {
        private static final long serialVersionUID = 1L;
        @Override public boolean test(final GATKRead read){return read.getStart() >= 1;}
    };
    static final ReadFilter endOk = new ReadFilter() {
        private static final long serialVersionUID = 1L;
        @Override public boolean test(final GATKRead read){return read.getEnd() <= 10;}
    };

    @DataProvider(name = "readsStartEnd")
    public Object[][] readsStartEnd(){
        return new Object[][]{
                { goodRead, true, true},
                { startBad, false, true},
                { endBad, true, false},
                { bothBad, false, false}
        };
    }


    @Test(dataProvider = "readsStartEnd")
    public void testTest(GATKRead read, boolean start, boolean end){
        Assert.assertEquals(startOk.test(read), start);
        Assert.assertEquals(endOk.test(read), end);
    }

    @Test(dataProvider = "readsStartEnd")
    public void testNegate(GATKRead read, boolean start, boolean end){
        Assert.assertEquals(startOk.negate().test(read), !start);
        Assert.assertEquals(endOk.negate().test(read), !end);
    }

    @DataProvider(name = "readsAnd")
    public Object[][] readsAnd(){
        return new Object[][]{
                { goodRead, true},
                { startBad, false},
                { endBad, false},
                { bothBad, false}
        };
    }

    @Test(dataProvider = "readsAnd")
    public void testAnd(GATKRead read, boolean expected){
        boolean aTest = startOk.test(read);
        boolean bTest = endOk.test(read);
        boolean cTest = endOk.negate().test(read);

        ReadFilter startAndEndOk = startOk.and(endOk);
        ReadFilter endAndStartOk = endOk.and(startOk);
        Assert.assertEquals(startAndEndOk.test(read), expected);
        Assert.assertEquals(endAndStartOk.test(read), expected);
    }

    @DataProvider(name = "readsOr")
    public Object[][] readsOr(){
        return new Object[][]{
                { goodRead, true},
                { startBad, true},
                { endBad, true},
                { bothBad, false}
        };
    }


    @Test(dataProvider = "readsOr")
    public void testOr(GATKRead read, boolean expected) {
        ReadFilter startOrEndOk = startOk.or(endOk);
        ReadFilter endOrStartOk = endOk.or(startOk);
        Assert.assertEquals(startOrEndOk.test(read), expected);
        Assert.assertEquals(endOrStartOk.test(read), expected);
    }

    @DataProvider(name = "deeper")
    public Object[][] deeper(){
        return new Object[][]{
                { goodRead, false},
                { startBad, true},
                { endBad, true},
                { bothBad, false}
        };
    }

    @Test(dataProvider = "deeper")
    public void testDeeperChaining(GATKRead read, boolean expected){
        ReadFilter notAMinionOfGozer = new ReadFilter() {
            private static final long serialVersionUID = 1L;
            @Override public boolean test(final GATKRead read){return !read.getName().equals("Zuul");}
        };
        ReadFilter readChecksOut = startOk.or(endOk).and(notAMinionOfGozer);
        Assert.assertEquals(readChecksOut.test(read), expected);
        Assert.assertEquals(readChecksOut.and(readChecksOut).test(read), expected);
        Assert.assertEquals(readChecksOut.and(r -> false).test(read), false);
        Assert.assertEquals(readChecksOut.or(r -> true).test(read), true);
    }

    @Test
    public void testFromListNull() {
        ReadFilter rf = ReadFilter.fromList(null, ArtificialReadUtils.createArtificialSamHeader(1, 1, 10));
        Assert.assertTrue(rf.getClass() == ReadFilterLibrary.AllowAllReadsReadFilter.class);
    }
    @Test
    public void testFromListEmpty() {
        ReadFilter rf = ReadFilter.fromList(Collections.emptyList(), ArtificialReadUtils.createArtificialSamHeader(1, 1, 10));
        Assert.assertTrue(rf.getClass() == ReadFilterLibrary.AllowAllReadsReadFilter.class);
    }

    @Test
    public void testFromListSingle() {
        List<ReadFilter> filters = new ArrayList<>();
        filters.add(ReadFilterLibrary.MAPPED);
        ReadFilter rf = ReadFilter.fromList(filters, ArtificialReadUtils.createArtificialSamHeader(1, 1, 10));
        Assert.assertTrue(rf.getClass() == ReadFilterLibrary.MAPPED.getClass());
    }

    @Test
    public void testFromListMultiOrdered() {
        List<ReadFilter> filters = new ArrayList<>();
        filters.add(ReadFilterLibrary.MAPPING_QUALITY_AVAILABLE);
        filters.add(ReadFilterLibrary.MAPPED);
        filters.add(ReadFilterLibrary.GOOD_CIGAR);
        filters.add(ReadFilterLibrary.FIRST_OF_PAIR);

        // Since we want to ensure that order of the input is honored, we need to test the
        // structure of the filter rather than the result
        ReadFilter rf = ReadFilter.fromList(filters, ArtificialReadUtils.createArtificialSamHeader(1, 1, 10));

        int count = verifyAndFilterOrder(
                rf,
                new String[] {
                        ReadFilterLibrary.MAPPING_QUALITY_AVAILABLE.getClass().getSimpleName(),
                        ReadFilterLibrary.MAPPED.getClass().getSimpleName(),
                        ReadFilterLibrary.GOOD_CIGAR.getClass().getSimpleName(),
                        ReadFilterLibrary.FIRST_OF_PAIR.getClass().getSimpleName()
                }
        );
        Assert.assertEquals(count, 4);
    }

    public static int verifyAndFilterOrder(final ReadFilter rf, final String[] expectedOrder) {
        Assert.assertEquals(rf.getClass(), ReadFilter.ReadFilterAnd.class);
        return verifyAndFilterOrder((ReadFilter.ReadFilterAnd) rf, expectedOrder);
    }

    public static int verifyAndFilterOrder(final ReadFilter.ReadFilterAnd rf, final String[] expectedOrder) {
        if (rf.lhs instanceof ReadFilter.ReadFilterAnd) {
            int count = verifyAndFilterOrder((ReadFilter.ReadFilterAnd) rf.lhs, expectedOrder);
            Assert.assertEquals(expectedOrder[count], rf.rhs.getClass().getSimpleName());
            return ++count;
        } else {
            int count = 0;
            Assert.assertEquals(expectedOrder[count], rf.lhs.getClass().getSimpleName());
            count++;
            Assert.assertEquals(expectedOrder[count], rf.rhs.getClass().getSimpleName());
            count++;
            return count;
        }
    }

}