package org.broadinstitute.hellbender.tools.picard.illumina.parser.readers;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;

public class AbstractIlluminaPositionFileReaderTest {

    //Pos values should never be negative
    private static final float [] X_COORDS = new float[] {
        0.00f,      0.012f, 0.09f,   1.44f,       23.45f,           155.893f,  588.921f, 801f,    1201.182f, 101.25f,
        13201.01f,  0.59f,  0.3540f, 9999999.99f, 7989.88999f,      9298.134f, 12.3f,    109.09f, 54.3f,     17.15f
    };

    private static final int [] X_QSEQ_COORDS = new int [] {
        1000,   1000, 1001, 1014,      1235,       2559,  6889, 9010, 13012, 2013,
        133010, 1006, 1004, 100001000, 80899,      93981, 1123, 2091, 1543, 1172
    };

    private static final float [] Y_COORDS = new float[] {
        //A note on index 2 -> this number causes truncation of the float, presumably that's fine as
        //Illumina describes the value as a floating point number
        679.99f,    32.0145f, 1878854.89f, 0.00f, 23.45f,           9875.64f,  42.42f,  64.01530f, 103.25f, 102.75f,
        13201.01f,  0.59f,    0.3540f,     0.02f, 9875.66f,         9298.134f, 12.3f,   67.012f,   0.1254f, 19.54f
    };

    private static final int [] Y_QSEQ_COORDS = new int [] {
        7800,   1320, 18789548, 1000, 1235,       99756, 1424, 1640, 2033, 2028,
        133010, 1006, 1004,     1000, 99757,      93981, 1123, 1670, 1001, 1195
    };


    private static class MockPositionFileReader extends AbstractIlluminaPositionFileReader {
        private final float [] xCoords;
        private final float [] yCoords;
        private int currentCluster;

        public MockPositionFileReader(final String fileName, final int lane, final int tile, final float [] xCoords, final float [] yCoords) {
            super(new File(fileName));
            this.xCoords = xCoords;
            this.yCoords = yCoords;
            currentCluster = 0;
        }

        public MockPositionFileReader(final int lane, final int tile, final float [] xCoords, final float [] yCoords) {
            super(new File("s_" + lane + "_" + tile + "_pos.txt.gz"));
            this.xCoords = xCoords;
            this.yCoords = yCoords;
            currentCluster = 0;
        }

        @Override
        protected PositionInfo unsafeNextInfo() {
            PositionInfo pi = new PositionInfo(xCoords[currentCluster], yCoords[currentCluster], getLane(), getTile());
            ++currentCluster;
            return pi;
        }

        @Override
        protected String makeExceptionMsg() {
            return "Abstract IlluminaPositionFileReaderTest currentCluster=" + currentCluster;
        }

        @Override
        public boolean hasNext() {
            return currentCluster < xCoords.length;
        }

        @Override
        public void close() {
        }
    }

    @Test
    public void validReaderTest() {
        final int lane = 2;
        final int tile = 8;
        final AbstractIlluminaPositionFileReader reader = new MockPositionFileReader(lane, tile, X_COORDS, Y_COORDS);

        int index = 0;
        while(reader.hasNext()) {
            AbstractIlluminaPositionFileReader.PositionInfo pi = reader.next();
            Assert.assertEquals(pi.lane, lane);
            Assert.assertEquals(pi.tile, tile);
            Assert.assertEquals(pi.xPos, X_COORDS[index]);
            Assert.assertEquals(pi.yPos, Y_COORDS[index]);
            Assert.assertEquals(pi.xQseqCoord, X_QSEQ_COORDS[index]);
            Assert.assertEquals(pi.yQseqCoord, Y_QSEQ_COORDS[index]);
            ++index;
        }

        if(index < X_COORDS.length) {
            throw new RuntimeException("Elements still remaining in test data!");
        }
    }

    @DataProvider(name="invalidPositions")
    public Object[][] invalidPositions() {
        return new Object[][] {
            {new float[]{0f, 5f, -11f},    new float[]{1f, 12f, 13f}},
            {new float[]{-15.05f},         new float[]{19801f}},
            {new float[]{10.05f, 3f, 8f}, new float[]{-120899.723f, 4f, 9f}},
            {new float[]{9.0f, 2.3f, AbstractIlluminaPositionFileReader.MAX_POS + 1}, new float[]{3.2f, 8.1f, 99.1f}},
            {new float[]{0.01f}, new float[]{AbstractIlluminaPositionFileReader.MAX_POS + 1000}}
        };
    }

    @Test(expectedExceptions = IllegalArgumentException.class, dataProvider = "invalidPositions")
    public void invalidReaderTest(float [] xCoords, float [] yCoords){
        final int lane = 3;
        final int tile = 4;

        final AbstractIlluminaPositionFileReader reader = new MockPositionFileReader(lane, tile, xCoords, yCoords);
        int index = 0;
        while(reader.hasNext()) {
            AbstractIlluminaPositionFileReader.PositionInfo pi = reader.next();
            Assert.assertEquals(pi.lane, lane);
            Assert.assertEquals(pi.tile, tile);
            Assert.assertEquals(pi.xPos, xCoords[index]);
            Assert.assertEquals(pi.yPos, yCoords[index]);
            ++index;
        }
    }

    @DataProvider(name = "invalidFileNames")
    public Object[][] invalidFileNames() {
        return new Object[][]{
            {"whatever.locs"},
            {"whatever.clocs"},
            {"whatever.pos"},
            {"s_1.clocs"},
            {"s__2.clocs"},
            {"s_1_4.Notclocs"},
            {"a_1_4.pos"},
            {"a_1_4.pos.txt"}
        };
    }

    @Test(expectedExceptions = IlluminaReaderException.class, dataProvider = "invalidFileNames")
    public void invalidFileNamesTest(final String fileName){
        final int lane = 3;
        final int tile = 4;

        final AbstractIlluminaPositionFileReader reader = new MockPositionFileReader(fileName, 0, 0, null, null);
    }

    @Test(expectedExceptions = UnsupportedOperationException.class)
    public void iteratorRemoveTest() {
        final AbstractIlluminaPositionFileReader reader = new MockPositionFileReader("s_1_1_pos.txt", 0, 0, null, null);
        reader.remove();
    }
}
