package org.broadinstitute.hellbender.tools.picard.illumina.parser.readers;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;

public final class LocsFileReaderTest {
    private static final File TestDir = new File("src/test/resources/org/broadinstitute/hellbender/tools/picard/illumina/readerTests");
    public  static final File LocsFile = new File(TestDir, "s_1_6.locs");
    public  static final int ExpectedTile = 6;
    public  static final int ExpectedLane = 1;
    public  static final int NumValues = 200;
    public  static final float [][] FloatCoords = {
        {1703.0117f, 64.01593f}, {1660.3038f, 64.08882f},  {1769.7501f, 65.12467f}, {1726.6725f, 68.367805f}, {1401.213f,  72.07282f},
        {1358.2775f, 72.07892f}, {1370.5197f, 77.699715f}, {1661.5403f, 77.70719f}, {1682.0504f, 78.76725f},  {1563.8765f, 79.08009f}
    };
    public static final int [][] QSeqCoords = {
        {18030, 1640}, {17603, 1641}, {18698, 1651}, {18267, 1684}, {15012, 1721},
        {14583, 1721}, {14705, 1777}, {17615, 1777}, {17821, 1788}, {16639, 1791}
    };
    public static final int [] Indices = {
        0,     1,  19,  59,  100,
        101,  179, 180, 198, 199
    };

    @Test
    public void passingFileTest() {
        final LocsFileReader reader = new LocsFileReader(LocsFile);

        int tdIndex = 0;
        int nextIndex = Indices[tdIndex];
        for(int i = 0; i < NumValues; i++) {
            Assert.assertTrue(reader.hasNext());

            final AbstractIlluminaPositionFileReader.PositionInfo piLocs = reader.next();

            if(i == nextIndex) {
                PosFileReaderTest.comparePositionInfo(piLocs, FloatCoords[tdIndex][0], FloatCoords[tdIndex][1],
                                                              QSeqCoords[tdIndex][0],  QSeqCoords[tdIndex][1],
                                                              ExpectedLane, ExpectedTile, i);

                if(tdIndex < Indices.length-1) {
                    nextIndex = Indices[++tdIndex];
                }
            }
        }

        Assert.assertFalse(reader.hasNext());
    }

    @DataProvider(name = "invalidFiles")
    public Object[][]invalidFiles() {
        return new Object[][] {
            {"s_1_7.locs"},
            {"s_1_8.locs"},
            {"s_1_9.locs"},
            {"s_1_10.locs"},
            {"s_f2af.locs"}
        };
    }

    @Test(expectedExceptions = IlluminaReaderException.class, dataProvider = "invalidFiles")
    public void invalidFilesTest(final String fileName) {
        final LocsFileReader reader = new LocsFileReader(new File(TestDir, fileName));
    }
}
