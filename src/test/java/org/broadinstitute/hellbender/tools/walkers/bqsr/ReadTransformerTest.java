package org.broadinstitute.hellbender.tools.walkers.bqsr;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordFactory;
import htsjdk.samtools.SAMUtils;
import org.broadinstitute.hellbender.utils.sam.ArtificialSAMUtils;
import org.testng.Assert;
import org.testng.annotations.BeforeTest;
import org.testng.annotations.Test;

import static org.testng.Assert.*;

public class ReadTransformerTest {
    static SAMFileHeader header =ArtificialSAMUtils.createArtificialSamHeader();

    SAMRecord sam;

    static ReadTransformer moveLeft = samrecord ->  {
        samrecord.setAlignmentStart(samrecord.getAlignmentStart() - 1);
        return samrecord;
        };

    static ReadTransformer moveTo20 = samrecord -> {
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
    }

    @Test
    public void testCompose() throws Exception {
        SAMRecord sam = ArtificialSAMUtils.createArtificialRead(header, "read1",0, 10,10);
        ReadTransformer leftThen20 = moveLeft.compose(moveTo20);
        leftThen20.apply(sam);
        Assert.assertEquals(sam.getAlignmentStart(), 19);
    }

    @Test
    public void testIdentity() throws Exception {
        SAMRecord sam = ArtificialSAMUtils.createArtificialRead(header, "read1",0, 10,10);
        Assert.assertEquals(ReadTransformer.identity().apply(sam), sam);

    }
}