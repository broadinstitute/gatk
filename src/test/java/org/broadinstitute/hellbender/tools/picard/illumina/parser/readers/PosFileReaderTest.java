package org.broadinstitute.hellbender.tools.picard.illumina.parser.readers;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;

public final class PosFileReaderTest {
    //need to test for negatives
    public static final File TestDir = new File("src/test/resources/org/broadinstitute/hellbender/tools/picard/illumina/readerTests");
    public static final File InvalidNamePosFile = new File(TestDir, "pos_failing1_pos.txt");

    public static final File PassingPosFile = new File(TestDir, "s_2_1101_pos.txt");
    public static final int PassingTile = 1101;
    public static final int PassingLane = 2;
    public static final float[][] PassingPosFloatCoord = {
            {17.80f, 0.30f}, {23.70f, 1.70f}, {18.50f, 3.90f}, {22.80f, 4.00f}, {15.90f, 4.20f},
            {19.10f, 5.60f}, {23.90f, 6.70f}, {16.70f, 7.00f}, {21.50f, 7.10f}, {18.90f, 9.30f},
            {17.10f, 10.30f}, {19.90f, 10.80f}, {15.90f, 11.90f}, {21.60f, 12.00f}, {24.80f, 12.40f},
            {17.90f, 14.20f}, {22.00f, 14.50f}, {23.60f, 15.20f}, {20.00f, 15.50f}, {16.30f, 17.30f}
    };
    public static final int[][] PassingPosQSeqCoord = {
            {1178, 1003}, {1237, 1017}, {1185, 1039}, {1228, 1040}, {1159, 1042},
            {1191, 1056}, {1239, 1067}, {1167, 1070}, {1215, 1071}, {1189, 1093},
            {1171, 1103}, {1199, 1108}, {1159, 1119}, {1216, 1120}, {1248, 1124},
            {1179, 1142}, {1220, 1145}, {1236, 1152}, {1200, 1155}, {1163, 1173}
    };

    @Test
    public void validPosFileTest() {
        final PosFileReader pfr = new PosFileReader(PassingPosFile);
        for (int i = 0; i < PassingPosFloatCoord.length; i++) {
            Assert.assertTrue(pfr.hasNext());
            final AbstractIlluminaPositionFileReader.PositionInfo pi = pfr.next();
            comparePositionInfo(pi, i);
        }
        Assert.assertFalse(pfr.hasNext());
        pfr.close();
    }

    public static void comparePositionInfo(final AbstractIlluminaPositionFileReader.PositionInfo pi, final int index) {
        comparePositionInfo(pi, PassingPosFloatCoord[index][0], PassingPosFloatCoord[index][1],
                PassingPosQSeqCoord[index][0], PassingPosQSeqCoord[index][1],
                PassingLane, PassingTile, index);
    }

    public static void comparePositionInfo(final AbstractIlluminaPositionFileReader.PositionInfo pi, final float xPos, final float yPos,
                                           final int xQSeqCoord, final int yQSeqCoord, final int lane, final int tile, final int index) {
        Assert.assertEquals(pi.xPos, xPos, "Differs at index: " + index);
        Assert.assertEquals(pi.yPos, yPos, "Differs at index: " + index);
        Assert.assertEquals(pi.xQseqCoord, xQSeqCoord, "Differs at index: " + index);
        Assert.assertEquals(pi.yQseqCoord, yQSeqCoord, "Differs at index: " + index);
        Assert.assertEquals(pi.lane, lane, "Differs at index: " + index);
        Assert.assertEquals(pi.tile, tile, "Differs at index: " + index);
    }

    @Test(expectedExceptions = IlluminaReaderException.class)
    public void failingNamePosFileTest() {
        final PosFileReader pfr = new PosFileReader(InvalidNamePosFile);
    }

    @DataProvider(name = "invalidDataFiles")
    public Object[][] invalidDataFiles() {
        return new Object[][]{
                {new File(TestDir, "s_1_1101_pos.txt")},
                {new File(TestDir, "s_1_1102_pos.txt")},
                {new File(TestDir, "s_1_1103_pos.txt")}
        };
    }

    @Test(dataProvider = "invalidDataFiles", expectedExceptions = IlluminaReaderException.class)
    public void failingDataPosFileTest(final File file) {
        final PosFileReader pfr = new PosFileReader(file);
        try {
            while (pfr.hasNext()) {
                pfr.next();
            }
        } finally {
            pfr.close();
        }
    }

    @Test
    public void zeroLengthFileTest() {
        final PosFileReader pfr = new PosFileReader(new File(TestDir, "s_1_1104_pos.txt"));
        Assert.assertFalse(pfr.hasNext());
        pfr.close();
    }
}
