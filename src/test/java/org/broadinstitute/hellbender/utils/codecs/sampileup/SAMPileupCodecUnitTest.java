package org.broadinstitute.hellbender.utils.codecs.sampileup;

import org.testng.Assert;
import org.testng.annotations.Test;

public class SAMPileupCodecUnitTest {

    @Test
    public void testCanDecode() {
        final String EXTRA_CHAR = "1";
        SAMPileupCodec codec = new SAMPileupCodec();
        for(final String ext: SAMPileupCodec.FILE_EXTENSIONS) {
            Assert.assertTrue(codec.canDecode("filename." + ext));
            Assert.assertTrue(codec.canDecode("filename" + EXTRA_CHAR + "." + ext));
            Assert.assertFalse(codec.canDecode("filename." + ext + "1"));
            Assert.assertFalse(codec.canDecode("filename" + ext));
        }
    }
}