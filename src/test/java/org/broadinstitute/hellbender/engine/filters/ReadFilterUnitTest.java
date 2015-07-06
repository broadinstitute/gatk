package org.broadinstitute.hellbender.engine.filters;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

public final class ReadFilterUnitTest {

    // Mirrors CountedReadFilterUnitTest
    static final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader(1, 1, 10);
    static final GATKRead goodRead = ArtificialReadUtils.createArtificialRead(header, "Zuul", 0, 2,2);
    static final GATKRead endBad = ArtificialReadUtils.createArtificialRead(header, "Peter", 0, 1,100);
    static final GATKRead startBad = ArtificialReadUtils.createArtificialRead(header, "Ray", 0, -1,2);
    static final GATKRead bothBad = ArtificialReadUtils.createArtificialRead(header, "Egon", 0, -1,100);
    static final ReadFilter startOk = r -> r.getStart() >= 1;
    static final ReadFilter endOk = r -> r.getEnd() <= 10;

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
        ReadFilter notAMinionOfGozer = r -> !r.getName().equals("Zuul");
        ReadFilter readChecksOut = startOk.or(endOk).and(notAMinionOfGozer);
        Assert.assertEquals(readChecksOut.test(read), expected);
        Assert.assertEquals(readChecksOut.and(readChecksOut).test(read), expected);
        Assert.assertEquals(readChecksOut.and(r -> false).test(read), false);
        Assert.assertEquals(readChecksOut.or(r -> true).test(read), true);
    }
}