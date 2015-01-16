package org.broadinstitute.hellbender.tools.walkers.bqsr;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.utils.sam.ArtificialSAMUtils;
import org.testng.Assert;
import org.testng.annotations.BeforeTest;
import org.testng.annotations.Test;

public class ReadTransformerTest {
    static final SAMFileHeader header =ArtificialSAMUtils.createArtificialSamHeader();
    SAMRecord sam;

    static final ReadTransformer moveLeft = samrecord ->  {
        samrecord.setAlignmentStart(samrecord.getAlignmentStart() - 1);
        return samrecord;
        };

    static final ReadTransformer moveTo20 = samrecord -> {
        samrecord.setAlignmentStart(20);
        return samrecord;
    };

    @BeforeTest
    private void resetSamRecord(){
        sam = ArtificialSAMUtils.createArtificialRead(header, "read1",0, 10,10);
    }

    @Test
    public void testAndThen() {
        ReadTransformer leftThen20 = moveLeft.andThen(moveTo20);
        leftThen20.apply(sam);
        Assert.assertEquals(sam.getAlignmentStart(), 20);
        moveTo20.andThen(moveLeft).apply(sam);
        Assert.assertEquals(sam.getAlignmentStart(), 19);
    }

    @Test
    public void testCompose(){
        ReadTransformer leftThen20 = moveLeft.compose(moveTo20);
        leftThen20.apply(sam);
        Assert.assertEquals(sam.getAlignmentStart(), 19);
        moveTo20.compose(moveLeft).apply(sam);
        Assert.assertEquals(sam.getAlignmentStart(), 20);
    }

    @Test
    public void testAndChain() {
        moveTo20.andThen(moveLeft)
                .andThen(moveLeft)
                .andThen(moveLeft)
                .andThen(moveLeft)
                .apply(sam);
        Assert.assertEquals(sam.getAlignmentStart(), 16 );
    }

    @Test
    public void testMixedChain(){

        moveLeft.compose(r -> { r.setMappingQuality(42); return r;})
                .compose(moveTo20)
                .andThen(moveLeft)
                .apply(sam);
        Assert.assertEquals(sam.getMappingQuality(), 42);
        Assert.assertEquals(sam.getAlignmentStart(), 18);

    }

    @Test
    public void testIdentity(){
        SAMRecord sam = ArtificialSAMUtils.createArtificialRead(header, "read1",0, 10,10);
        Assert.assertEquals(ReadTransformer.identity().apply(sam), sam);

    }
}