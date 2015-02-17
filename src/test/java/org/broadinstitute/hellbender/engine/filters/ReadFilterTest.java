package org.broadinstitute.hellbender.engine.filters;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.utils.sam.ArtificialSAMUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import static org.testng.Assert.*;

public class ReadFilterTest {

    static final SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(1, 1, 10);
    static final SAMRecord goodRead = ArtificialSAMUtils.createArtificialRead(header, "read1", 0, 2,2);
    static final SAMRecord badRead = ArtificialSAMUtils.createArtificialRead(header, "read2", 0, 1,100);
    static final ReadFilter startOk = r -> r.getAlignmentStart() >= 1;
    static final ReadFilter endOK = r -> r.getAlignmentEnd() <= 10;

    @Test
    public void testAnd(){
        ReadFilter startAndEndOk = startOk.and(endOK);
        Assert.assertTrue(startAndEndOk.test(goodRead));
        Assert.assertFalse(startAndEndOk.test(badRead));
    }

    @Test
    public void testOr() {
        ReadFilter startOrEndOk = startOk.or(endOK);
        Assert.assertTrue(startOrEndOk.test(goodRead));
        Assert.assertTrue(startOrEndOk.test(badRead));
    }

    @Test
    public void testNegate(){
        Assert.assertTrue(startOk.test(goodRead));
        Assert.assertFalse(startOk.negate().test(goodRead));
    }
}