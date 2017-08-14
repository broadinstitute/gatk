package org.broadinstitute.hellbender.utils.codecs;

import htsjdk.tribble.readers.LineIterator;
import htsjdk.tribble.readers.PositionalBufferedStream;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.test.BaseTest;

import org.broadinstitute.hellbender.utils.tsv.TableUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.ByteArrayInputStream;
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

    @DataProvider(name = "targetCodecMissingHeaders")
    public Object[][] targetCodecMissingHeaders() {
        return new Object[][] {
                { "" },
                { "# a comment only" },
                { "# a comment\n" +
                        "1" + TableUtils.COLUMN_SEPARATOR_STRING +
                        "389" + TableUtils.COLUMN_SEPARATOR_STRING +
                        "455" + TableUtils.COLUMN_SEPARATOR_STRING +
                        "target_1_389_455\n" },

        };
    }

    @Test(dataProvider = "targetCodecMissingHeaders", expectedExceptions = UserException.BadInput.class)
    public void testIndexableMissingHeaderFails(final String featureTableContents) {
        final TargetCodec codec = new TargetCodec();
        LineIterator lit = (LineIterator) codec.makeIndexableSourceFromStream(new ByteArrayInputStream(featureTableContents.getBytes()));
        codec.readActualHeader(lit);  // no header to be found
    }

    @DataProvider(name = "targetCodecFeatureOffsets")
    public Object[][] goodTargetCodecHeaders() {
        return new Object[][] {
                new Object[] {
                        "contig" + TableUtils.COLUMN_SEPARATOR_STRING +
                        "start" + TableUtils.COLUMN_SEPARATOR_STRING +
                        "stop" + TableUtils.COLUMN_SEPARATOR_STRING +
                        "name\n" +
                        "1" + TableUtils.COLUMN_SEPARATOR_STRING +
                        "389" + TableUtils.COLUMN_SEPARATOR_STRING +
                        "455" + TableUtils.COLUMN_SEPARATOR_STRING +
                        "target_1_389_455\n",
                        23L
                },
                new Object[] {
                        "# a comment\n" +
                        "contig" + TableUtils.COLUMN_SEPARATOR_STRING +
                        "start" + TableUtils.COLUMN_SEPARATOR_STRING +
                        "stop" + TableUtils.COLUMN_SEPARATOR_STRING +
                        "name\n" +
                        "1" + TableUtils.COLUMN_SEPARATOR_STRING +
                        "389" + TableUtils.COLUMN_SEPARATOR_STRING +
                        "455" + TableUtils.COLUMN_SEPARATOR_STRING +
                        "target_1_389_455\n",
                        35L
                },
                new Object[] {
                        "# one comment\n" +
                        "# two comments\n" +
                        "contig" + TableUtils.COLUMN_SEPARATOR_STRING +
                        "start" + TableUtils.COLUMN_SEPARATOR_STRING +
                        "stop" + TableUtils.COLUMN_SEPARATOR_STRING +
                        "name\n" +
                        "1" + TableUtils.COLUMN_SEPARATOR_STRING +
                        "389" + TableUtils.COLUMN_SEPARATOR_STRING +
                        "455" + TableUtils.COLUMN_SEPARATOR_STRING +
                        "target_1_389_455\n",
                        52L
                }
        };
    }

    @Test(dataProvider = "targetCodecFeatureOffsets")
    public void testIndexableFeatureOffset(final String featureTableContents, final long firstFeatureOffset) {
        final TargetCodec codec = new TargetCodec();
        final PositionalBufferedStream pbs = new PositionalBufferedStream(new ByteArrayInputStream(featureTableContents.getBytes()));
        final LineIterator lit = (LineIterator) codec.makeIndexableSourceFromStream(pbs);
        codec.readActualHeader(lit);
        Assert.assertEquals(pbs.getPosition(), firstFeatureOffset);
    }

}
