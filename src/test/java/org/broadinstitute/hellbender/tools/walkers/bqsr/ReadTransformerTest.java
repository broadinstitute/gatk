package org.broadinstitute.hellbender.tools.walkers.bqsr;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.transformers.ReadTransformer;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

public final class ReadTransformerTest {
    static final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader();
    
    static final ReadTransformer moveLeft = samrecord ->  {
        samrecord.setPosition(samrecord.getContig(), samrecord.getStart() - 1);
        return samrecord;
    };

    static final ReadTransformer moveTo20 = samrecord -> {
        samrecord.setPosition(samrecord.getContig(), 20);
        return samrecord;
    };

    @DataProvider(name = "UnmodifiedReadDataProvider")
    public Object[][] getUnmodifiedReadData() {
        return new Object[][] {
                { ArtificialReadUtils.createArtificialRead(header, "read1", 0, 10, 10) }
        };
    }

    @Test(dataProvider = "UnmodifiedReadDataProvider")
    public void testAndThen( final GATKRead read ) {
        ReadTransformer leftThen20 = moveLeft.andThen(moveTo20);
        leftThen20.apply(read);
        Assert.assertEquals(read.getStart(), 20);
        moveTo20.andThen(moveLeft).apply(read);
        Assert.assertEquals(read.getStart(), 19);
    }

    @Test(dataProvider = "UnmodifiedReadDataProvider")
    public void testCompose( final GATKRead read ) {
        ReadTransformer leftThen20 = moveLeft.compose(moveTo20);
        leftThen20.apply(read);
        Assert.assertEquals(read.getStart(), 19);
        moveTo20.compose(moveLeft).apply(read);
        Assert.assertEquals(read.getStart(), 20);
    }

    @Test(dataProvider = "UnmodifiedReadDataProvider")
    public void testAndChain( final GATKRead read ) {
        moveTo20.andThen(moveLeft)
                .andThen(moveLeft)
                .andThen(moveLeft)
                .andThen(moveLeft)
                .apply(read);
        Assert.assertEquals(read.getStart(), 16);
    }

    @Test(dataProvider = "UnmodifiedReadDataProvider")
    public void testMixedChain( final GATKRead read ) {

        moveLeft.compose(r -> { r.setMappingQuality(42); return r;})
                .compose(moveTo20)
                .andThen(moveLeft)
                .apply(read);
        Assert.assertEquals(read.getMappingQuality(), 42);
        Assert.assertEquals(read.getStart(), 18);

    }

    @Test(dataProvider = "UnmodifiedReadDataProvider")
    public void testIdentity( final GATKRead read ) {
        Assert.assertEquals(ReadTransformer.identity().apply(read), read);
    }
}