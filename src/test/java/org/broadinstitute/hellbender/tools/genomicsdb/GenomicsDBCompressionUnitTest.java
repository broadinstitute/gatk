package org.broadinstitute.hellbender.tools.genomicsdb;

import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.tools.genomicsdb.GenomicsDBCompression.Codec;
import org.testng.Assert;
import org.testng.annotations.Test;

public class GenomicsDBCompressionUnitTest extends GATKBaseTest {
    private static final String ARG = GenomicsDBImport.COMPRESSION_LONG_NAME;

    @Test
    public void codecWithoutLevelUsesLevelOne() {
        Assert.assertEquals(GenomicsDBCompression.parse(ARG, "zstd"), new GenomicsDBCompression(Codec.ZSTD, 1));
        Assert.assertEquals(GenomicsDBCompression.parse(ARG, "gzip"), new GenomicsDBCompression(Codec.GZIP, 1));
        Assert.assertEquals(GenomicsDBCompression.parse(ARG, "lz4"), new GenomicsDBCompression(Codec.LZ4, 1));
    }

    @Test
    public void noneWithoutLevelUsesLevelZero() {
        Assert.assertEquals(GenomicsDBCompression.parse(ARG, "none"), new GenomicsDBCompression(Codec.NONE, 0));
    }

    @Test
    public void codecWithLevelUsesThatLevel() {
        Assert.assertEquals(GenomicsDBCompression.parse(ARG, "zstd:7"), new GenomicsDBCompression(Codec.ZSTD, 7));
        Assert.assertEquals(GenomicsDBCompression.parse(ARG, "gzip:9"), new GenomicsDBCompression(Codec.GZIP, 9));
    }

    @Test
    public void codecNameIsCaseInsensitive() {
        Assert.assertEquals(GenomicsDBCompression.parse(ARG, "ZStd:2"), new GenomicsDBCompression(Codec.ZSTD, 2));
    }

    @Test
    public void codecsMapToTileDbCodes() {
        Assert.assertEquals(Codec.NONE.tiledbCode(), 0);
        Assert.assertEquals(Codec.GZIP.tiledbCode(), 1);
        Assert.assertEquals(Codec.ZSTD.tiledbCode(), 2);
        Assert.assertEquals(Codec.LZ4.tiledbCode(), 3);
    }

    @Test(expectedExceptions = CommandLineException.BadArgumentValue.class)
    public void unknownCodecIsRejected() {
        GenomicsDBCompression.parse(ARG, "brotli");
    }

    @Test(expectedExceptions = CommandLineException.BadArgumentValue.class)
    public void emptySpecificationIsRejected() {
        GenomicsDBCompression.parse(ARG, "");
    }

    @Test(expectedExceptions = CommandLineException.BadArgumentValue.class)
    public void nonIntegerLevelIsRejected() {
        GenomicsDBCompression.parse(ARG, "zstd:fast");
    }

    @Test(expectedExceptions = CommandLineException.BadArgumentValue.class)
    public void emptyLevelIsRejected() {
        GenomicsDBCompression.parse(ARG, "zstd:");
    }

    @Test(expectedExceptions = CommandLineException.BadArgumentValue.class)
    public void gzipLevelAboveNineIsRejected() {
        GenomicsDBCompression.parse(ARG, "gzip:10");
    }

    @Test(expectedExceptions = CommandLineException.BadArgumentValue.class)
    public void zstdLevelZeroIsRejected() {
        GenomicsDBCompression.parse(ARG, "zstd:0");
    }

    @Test(expectedExceptions = CommandLineException.BadArgumentValue.class)
    public void zstdLevelAboveTwentyTwoIsRejected() {
        GenomicsDBCompression.parse(ARG, "zstd:23");
    }

    @Test(expectedExceptions = CommandLineException.BadArgumentValue.class)
    public void levelForNoneIsRejected() {
        GenomicsDBCompression.parse(ARG, "none:1");
    }

    @Test(expectedExceptions = CommandLineException.BadArgumentValue.class)
    public void extraFieldIsRejected() {
        GenomicsDBCompression.parse(ARG, "zstd:1:2");
    }

    @Test
    public void whitespaceAroundCodecAndLevelIsIgnored() {
        Assert.assertEquals(GenomicsDBCompression.parse(ARG, " zstd : 3 "), new GenomicsDBCompression(Codec.ZSTD, 3));
    }

    @Test
    public void highestLevelOfEachCodecIsAccepted() {
        Assert.assertEquals(GenomicsDBCompression.parse(ARG, "gzip:9"), new GenomicsDBCompression(Codec.GZIP, 9));
        Assert.assertEquals(GenomicsDBCompression.parse(ARG, "zstd:22"), new GenomicsDBCompression(Codec.ZSTD, 22));
        Assert.assertEquals(GenomicsDBCompression.parse(ARG, "lz4:127"), new GenomicsDBCompression(Codec.LZ4, 127));
    }

    @Test(expectedExceptions = CommandLineException.BadArgumentValue.class)
    public void gzipLevelZeroIsRejected() {
        GenomicsDBCompression.parse(ARG, "gzip:0");
    }

    @Test(expectedExceptions = CommandLineException.BadArgumentValue.class)
    public void lz4LevelZeroIsRejected() {
        GenomicsDBCompression.parse(ARG, "lz4:0");
    }

    // TileDB stores a level as one signed byte, so an lz4 acceleration above 127 would silently wrap
    @Test(expectedExceptions = CommandLineException.BadArgumentValue.class)
    public void lz4LevelAbove127IsRejected() {
        GenomicsDBCompression.parse(ARG, "lz4:128");
    }

    @Test
    public void levelTooLargeForAnIntegerIsReportedWithTheRange() {
        try {
            GenomicsDBCompression.parse(ARG, "zstd:99999999999");
            Assert.fail("expected the level to be rejected");
        } catch (final CommandLineException.BadArgumentValue e) {
            Assert.assertTrue(e.getMessage().contains("integer between 1 and 22"), e.getMessage());
        }
    }
}
