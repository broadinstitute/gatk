package org.broadinstitute.hellbender.utils.codecs.refseq;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

public class RefSeqCodecUnitTest extends GATKBaseTest {

    @Test
    public void testCanDecode() {
        final String EXTRA_CHAR = "1";
        RefSeqCodec codec = new RefSeqCodec();
        Assert.assertTrue(codec.canDecode("filename." + RefSeqCodec.FILE_EXT));
        Assert.assertTrue(codec.canDecode("filename" + EXTRA_CHAR + "." + RefSeqCodec.FILE_EXT));
        Assert.assertFalse(codec.canDecode("filename." + RefSeqCodec.FILE_EXT + EXTRA_CHAR));
        Assert.assertFalse(codec.canDecode("filename" + RefSeqCodec.FILE_EXT));
    }
}
