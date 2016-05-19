package org.broadinstitute.hellbender.utils.codecs.sampileup;

import org.testng.Assert;
import org.testng.annotations.Test;

public class SAMPileupCodecUnitTest {

    @Test
    public void testCanDecode() {
        final String EXTRA_CHAR = "1";
        SAMPileupCodec codec = new SAMPileupCodec();
        Assert.assertTrue(codec.canDecode("filename." + SAMPileupCodec.FILE_EXT));
        Assert.assertTrue(codec.canDecode("filename" + EXTRA_CHAR + "." + SAMPileupCodec.FILE_EXT));
        Assert.assertFalse(codec.canDecode("filename." + SAMPileupCodec.FILE_EXT + "1"));
        Assert.assertFalse(codec.canDecode("filename" + SAMPileupCodec.FILE_EXT));
    }
}