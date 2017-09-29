package org.broadinstitute.hellbender.utils.codecs;

import org.broadinstitute.hellbender.utils.test.BaseTest;

import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;

public class TargetCodecUnitTest extends BaseTest {

    @Test
    public void testTargetCodecCanDecodeFileURIs() {
        final TargetCodec codec = new TargetCodec();
        final String tsvURI = new File(publicTestDir + "org/broadinstitute/hellbender/tools/exome/targets.tsv").toURI().toString();
        Assert.assertTrue(tsvURI.startsWith("file:"));
        Assert.assertTrue(codec.canDecode(tsvURI), "TargetCodec should be able to decode file URIs");
    }

    @Test
    public void testTargetCodecCanDecodeAbsoluteFilePath() {
        final TargetCodec codec = new TargetCodec();
        final String tsvPath = new File(publicTestDir + "org/broadinstitute/hellbender/tools/exome/targets.tsv").getAbsolutePath();
        Assert.assertTrue(codec.canDecode(tsvPath), "TargetCodec should be able to decode absolute file paths");
    }

    @Test
    public void testTargetCodecCanDecodeRelativeFilePath() {
        final TargetCodec codec = new TargetCodec();
        final String tsvPath = "src/test/resources/org/broadinstitute/hellbender/tools/exome/targets.tsv";
        Assert.assertTrue(codec.canDecode(tsvPath), "TargetCodec should be able to decode relative file paths");
    }
}
