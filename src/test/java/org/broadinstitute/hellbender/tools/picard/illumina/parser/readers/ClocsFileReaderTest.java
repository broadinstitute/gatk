package org.broadinstitute.hellbender.tools.picard.illumina.parser.readers;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;

public class ClocsFileReaderTest {

    private static final File TEST_DIR = new File("src/test/resources/org/broadinstitute/hellbender/tools/picard/illumina/readerTests");
    private static final File PASSING_CLOCS_FILE = new File(TEST_DIR, "s_2_1101.clocs");

    public static final int blockSize = 25;
    public static final File MULTI_BIN_PASSING_CLOCS_FILE = new File(TEST_DIR, "s_3_2102.clocs");
    public static final File MBCF_W_EMPTY_BINS_AT_START   = new File(TEST_DIR, "s_3_2103.clocs");
    public static final File MBCF_W_EMPTY_BINS_AT_END     = new File(TEST_DIR, "s_3_2104.clocs");
    public static final File MBCF_W_EMPTY_BINS_THROUGHOUT = new File(TEST_DIR, "s_3_2105.clocs");
    public static final int MULTI_BIN_PASSING_CLOCS_LANE = 3;
    public static final int MULTI_BIN_EXPECTED_NUM_VALUES = 271;

    public static final File MBCF_MULTI_ROW_FILE = new File(TEST_DIR, "s_1_2106.clocs");

    public static final float [][] MULTI_BIN_FLOAT_COORD = {
        {9.7f, 1f},      {16.9f, 1.3f}, {18.5f, 22.6f}, {39.1f, 0.5f}, {46.1f, 1f},
        {38.4f, 8.4f},   {59.7f, 1.1f}, {89.9f, 0.30f}, {87.1f, 0.6f}, {83.6f, 0.7f},
        {111.4f, 14.3f}, {109.1f, 15f}, {105.7f, 15f}
    };
    public static final int [][] MULTI_BIN_Q_SEQ_COORD = {
        {1097, 1010}, {1169, 1013}, {1185, 1226}, {1391, 1005}, {1461, 1010},
        {1384, 1084}, {1597, 1011}, {1899, 1003}, {1871, 1006}, {1836, 1007},
        {2114, 1143}, {2091, 1150}, {2057, 1150}
    };
    public static final int [] MULTI_BIN_INDICES = {
        0,     1,  34,  39,  40,
        59,  99, 163, 164, 165,
        268, 269, 270
    };


    public static final float [][] MULTI_BIN_FLOAT_COORD_MIXED_EMPTY_BINS = {
            {9.7f, 1f},      {16.9f, 1.3f},  {18.5f, 22.6f},  {89.1f, 0.5f},  {96.1f, 1f},
            {88.4f, 8.4f},   {109.7f, 1.1f}, {164.9f, 0.30f}, {162.1f, 0.6f}, {158.6f, 0.7f},
            {186.4f, 14.3f}, {184.1f, 15f},  {180.7f, 15f}
    };
    public static final int [][] MULTI_BIN_Q_SEQ_COORD_MIXED_EMPTY_BINS = {
            {1097, 1010}, {1169, 1013}, {1185, 1226}, {1891, 1005}, {1961, 1010},
            {1884, 1084}, {2097, 1011}, {2649, 1003}, {2621, 1006}, {2586, 1007},
            {2864, 1143}, {2841, 1150}, {2807, 1150}
    };

    //These are all essentially the same file but with 0 or more empty bins spread through them
    @DataProvider(name = "multiBinPassingClocsFiles")
    public Object[][] multiBinPassingClocsFiles() {
        return new Object[][] {
            {MULTI_BIN_PASSING_CLOCS_FILE, 0, 2102},
            {MBCF_W_EMPTY_BINS_AT_START,   2, 2103},
            {MBCF_W_EMPTY_BINS_AT_END,     0, 2104}
        };
    }

    @Test(dataProvider = "multiBinPassingClocsFiles")
    public void multiBinTest(final File multiBinPassingClocsFile, final int binShift, final int tile) {
        final ClocsFileReader clocsReader = new ClocsFileReader(multiBinPassingClocsFile);

        int tdIndex = 0;
        int nextIndex = MULTI_BIN_INDICES[tdIndex];
        for(int i = 0; i < MULTI_BIN_EXPECTED_NUM_VALUES; i++) {
            Assert.assertTrue(clocsReader.hasNext());
            final AbstractIlluminaPositionFileReader.PositionInfo piClocs = clocsReader.next();

            if(i == nextIndex) {
                PosFileReaderTest.comparePositionInfo(piClocs, MULTI_BIN_FLOAT_COORD[tdIndex][0] + binShift * blockSize,       MULTI_BIN_FLOAT_COORD[tdIndex][1],
                                                               MULTI_BIN_Q_SEQ_COORD[tdIndex][0] + binShift * blockSize * 10,  MULTI_BIN_Q_SEQ_COORD[tdIndex][1],
                                                               MULTI_BIN_PASSING_CLOCS_LANE, tile,
                                                               tdIndex);
                if(tdIndex < MULTI_BIN_INDICES.length-1) {
                    nextIndex = MULTI_BIN_INDICES[++tdIndex];
                }
            }
        }

        Assert.assertFalse(clocsReader.hasNext());
    }

    @Test
    public void multiBinMixedEmptyBinTest() {
        final ClocsFileReader clocsReader = new ClocsFileReader(MBCF_W_EMPTY_BINS_THROUGHOUT);

        int tdIndex = 0;
        int nextIndex = MULTI_BIN_INDICES[tdIndex];
        for(int i = 0; i < MULTI_BIN_EXPECTED_NUM_VALUES; i++) {
            Assert.assertTrue(clocsReader.hasNext());
            final AbstractIlluminaPositionFileReader.PositionInfo piClocs = clocsReader.next();

            if(i == nextIndex) {
                PosFileReaderTest.comparePositionInfo(piClocs, MULTI_BIN_FLOAT_COORD_MIXED_EMPTY_BINS[tdIndex][0],  MULTI_BIN_FLOAT_COORD_MIXED_EMPTY_BINS[tdIndex][1],
                                                               MULTI_BIN_Q_SEQ_COORD_MIXED_EMPTY_BINS[tdIndex][0],  MULTI_BIN_Q_SEQ_COORD_MIXED_EMPTY_BINS[tdIndex][1],
                                                               MULTI_BIN_PASSING_CLOCS_LANE, 2105, tdIndex);
                if(tdIndex < MULTI_BIN_INDICES.length-1) {
                    nextIndex = MULTI_BIN_INDICES[++tdIndex];
                }
            }
        }

        Assert.assertFalse(clocsReader.hasNext());
    }

    public static final int [] MULTI_ROW_INDICES = new int[]{
        1, 10, 200, 2000, 10000, 30000
    };
    public static final float [][] MULTI_ROW_FLOAT_COORD = new float[][] {
            {23.2f,  0.2f},    {16.7f,  5.7f},  {60.9f, 15.0f},
            {596.2f, 19.6f},   {999.3f, 47.3f}, {1015f, 115.5f}
    };

    public static final int [][] MULTI_ROW_Q_SEQ_COORD = new int[][] {
        {1232, 1002},  {1167, 1057 },  {1609, 1150},
        {6962, 1196},  {10993, 1473}, {11150, 2155}
    };

    @Test
    public void multiRowFileTest() {
        final ClocsFileReader clocsReader = new ClocsFileReader(MBCF_MULTI_ROW_FILE);

        int tdIndex = 0;
        int nextIndex = MULTI_ROW_INDICES[tdIndex];
        for(int i = 1; i <= 39891; i++) {
            Assert.assertTrue(clocsReader.hasNext(), " i == " + i);
            final AbstractIlluminaPositionFileReader.PositionInfo piClocs = clocsReader.next();

            if(i == nextIndex) {
                PosFileReaderTest.comparePositionInfo(piClocs, MULTI_ROW_FLOAT_COORD[tdIndex][0], MULTI_ROW_FLOAT_COORD[tdIndex][1],
                                                               MULTI_ROW_Q_SEQ_COORD[tdIndex][0], MULTI_ROW_Q_SEQ_COORD[tdIndex][1],
                                                               1, 2106, tdIndex);

                if(tdIndex < MULTI_ROW_INDICES.length-1) {
                    nextIndex = MULTI_ROW_INDICES[++tdIndex];
                }
            }
        }

        Assert.assertFalse(clocsReader.hasNext());

    }

    @Test
    public void singleBinTest() {
        final ClocsFileReader clocsReader = new ClocsFileReader(PASSING_CLOCS_FILE);

        for(int i = 0; i < PosFileReaderTest.PassingPosFloatCoord.length; i++) {
            Assert.assertTrue(clocsReader.hasNext());
            final AbstractIlluminaPositionFileReader.PositionInfo piClocs = clocsReader.next();
            PosFileReaderTest.comparePositionInfo(piClocs, i);
        }

        Assert.assertFalse(clocsReader.hasNext());
    }
}
