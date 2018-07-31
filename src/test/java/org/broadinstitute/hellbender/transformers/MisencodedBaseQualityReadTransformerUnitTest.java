package org.broadinstitute.hellbender.transformers;


import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;

/**
 * Basic unit test for misencoded quals
 */
public final class MisencodedBaseQualityReadTransformerUnitTest extends GATKBaseTest {

    private SAMFileHeader header;

    @BeforeMethod
    public void before() {
        header = ArtificialReadUtils.createArtificialSamHeader(1, 1, 1000);
    }

    private GATKRead createRead(final byte[] quals) {
        final String readBases = "AAAAAAAAAA";
        GATKRead read = ArtificialReadUtils.createArtificialRead(header, "foo", 0, 10, readBases.getBytes(), quals);
        read.setCigar("10M");
        return read;
    }

    @Test
    public void testGoodQuals() {
        final byte[] goodQuals = { 60, 60, 60, 60, 60, 60, 60, 60, 60, 60 };
        final ReadTransformer tr = new MisencodedBaseQualityReadTransformer();
        GATKRead read = createRead(goodQuals);
        GATKRead newRead = tr.apply(read);
        Assert.assertEquals(read, newRead);
    }

    @Test
    public void testFixBadQuals() {
        final byte[] fixedQuals = { 28, 29, 31, 32, 33, 30, 31, 27, 26, 25 };
        final byte[] badQuals = { 59, 60, 62, 63, 64, 61, 62, 58, 57, 56 };
        final ReadTransformer tr = new MisencodedBaseQualityReadTransformer();
        final GATKRead read = createRead(badQuals);
        final GATKRead fixedRead = tr.apply(read);
        Assert.assertEquals(fixedQuals, fixedRead.getBaseQualities());
    }

    @Test(expectedExceptions = UserException.BadInput.class)
    public void testFixGoodQualsBlowUp() {
        final byte[] fixedQuals = { 28, 29, 31, 32, 33, 30, 31, 27, 26, 25 };
        final ReadTransformer tr = new MisencodedBaseQualityReadTransformer();
        final GATKRead read = createRead(fixedQuals);
        tr.apply(read);
    }
}