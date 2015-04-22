package org.broadinstitute.hellbender.tools.picard.illumina.parser;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.*;

import static htsjdk.samtools.util.CollectionUtil.makeList;

public final class PosParserTest {

    public static final File TEST_DIR = new File("src/test/resources/org/broadinstitute/hellbender/tools/picard/illumina/parserTests/posParser");
    public static final File S_1_1_POS = new File(TEST_DIR, "s_1_1_pos.txt");
    public static final File S_1_2_POS = new File(TEST_DIR, "s_1_2_pos.txt");
    public static final File S_1_3_POS = new File(TEST_DIR, "s_1_3_pos.txt");
    public static final File S_9_3201_POS = new File(TEST_DIR, "s_9_3201_pos.txt");
    public static final File S_9_3202_POS = new File(TEST_DIR, "s_9_3202_pos.txt");

    static class TestResult {
        public final int lane;
        public final int tile;
        public final float xPos;
        public final float yPos;
        public final int xQSeqCoord;
        public final int yQSeqCoord;

        public TestResult(final int lane, final int tile, final float xPos, final float yPos, final int xQseqCoord, final int yQSeqCoord) {
            this.lane = lane;
            this.tile = tile;
            this.xPos = xPos;
            this.yPos = yPos;
            this.xQSeqCoord = xQseqCoord;
            this.yQSeqCoord = yQSeqCoord;
        }
    }

    public static List<TestResult> makeTestResults(final int lane, final int tile, final float [] xyPos, final int [] xyQseq) {
        final ArrayList<TestResult> results = new ArrayList<TestResult>();
        for(int i = 0; i < xyPos.length; i+=2) {
            results.add(new TestResult(lane, tile, xyPos[i], xyPos[i+1], xyQseq[i], xyQseq[i+1]));
        }
        return results;
    }

    public static Map<String, List<TestResult>> TEST_DATA = new HashMap<String, List<TestResult>>();
    static {
       float[] pos = {
               101.35f,      207.8f,
               102.88f,      209.22f,
               0.35f,        211.99f,
               10.11f,       2540.55f,
               11011.81f,    211
       };

       int[] qseq = {
               2014,      3078,
               2029,      3092,
               1004,      3120,
               1101,      26406,
               111118,    3110
       };
       TEST_DATA.put(S_1_1_POS.getName(), makeTestResults(1,1, pos, qseq));

        pos = new float[] {
                2101.35f,  207.82f,
                202.88f,  209.222f,
                2.35f,    211.992f,
                2.11f,   2540.552f,
                21011.81f,    2112f
        };

        qseq = new int[]{
                22014,     3078,
                3029,      3092,
                1024,      3120,
                1021,      26406,
                211118,    22120
        };
        TEST_DATA.put(S_1_2_POS.getName(), makeTestResults(1,2, pos, qseq));

        pos = new float[] {
                301.35f,  207.83f,
                302.88f,  209.23f,
                0.35f,    911.993f,
                30.11f,   2540.553f
        };

        qseq = new int[]{
                4014,      3078,
                4029,      3092,
                1004,      10120,
                1301,      26406
        };
        TEST_DATA.put(S_1_3_POS.getName(), makeTestResults(1,3, pos, qseq));

        pos = new float[] {
            901.00f,   8011.00f,
            9.05f,    32.00f,
            3201.11f, 32001.00f
        };

        qseq = new int[]{
            10010, 81110,
            1091,  1320,
            33011, 321010
        };
        TEST_DATA.put(S_9_3201_POS.getName(), makeTestResults(9,3201, pos, qseq));

        pos = new float[] {
                3252f,      7785f,
                97.350f,    0.01f,
                88.00f,     942.01f,
                87.85f,     70.01f,
                32022.000f, 32000.0f
        };

        qseq = new int[]{
                33520,  78850,
                1974,   1000,
                1880,   10420,
                1879,   1700,
                321220, 321000
        };
        TEST_DATA.put(S_9_3202_POS.getName(), makeTestResults(9,3202, pos, qseq));
    }

    public static void compareResults(final TestResult tr, final PositionalData pd, final int index) {
        Assert.assertEquals(tr.xQSeqCoord, pd.getXCoordinate(), " At index " + index);
        Assert.assertEquals(tr.yQSeqCoord, pd.getYCoordinate(), " At index " + index);
    }

    @DataProvider(name = "singleTileData")
    public Object[][] singleTileData() {
        return new Object[][]{
            {1, 1, S_1_1_POS},
            {1, null, S_1_1_POS},
            {3202, 3202, S_9_3202_POS},
            {3202, null, S_9_3202_POS}
        };
    }

    @Test(dataProvider = "singleTileData")
    public void singleTileDataTest(final int tile, final Integer startingTile, final File file) {
        final IlluminaFileMap fm = new IlluminaFileMap();
        fm.put(tile, file);

        final PosParser parser = (startingTile == null) ? new PosParser(fm, IlluminaFileUtil.SupportedIlluminaFormat.Pos) :
                                                          new PosParser(fm, startingTile, IlluminaFileUtil.SupportedIlluminaFormat.Pos);
        final Iterator<TestResult> expected = TEST_DATA.get(file.getName()).iterator();

        int index = 0;
        while(expected.hasNext()) {
            Assert.assertTrue(parser.hasNext());
            compareResults(expected.next(), parser.next(), index);
            ++index;
        }

        Assert.assertFalse(parser.hasNext());
        parser.close();
    }

    @DataProvider(name = "multiTileData")
    public Object[][] multiTileData() {
        return new Object[][]{
            {makeList(1,2,3),     null,    makeList(S_1_1_POS, S_1_2_POS, S_1_3_POS)},
            {makeList(1,2,3),     1,       makeList(S_1_1_POS, S_1_2_POS, S_1_3_POS)},
            {makeList(1,3),       0,       makeList(S_1_1_POS, S_1_3_POS)},
            {makeList(3201,3202), null,    makeList(S_9_3201_POS, S_9_3202_POS)},
            {makeList(3201,3202), 1,       makeList(S_9_3201_POS, S_9_3202_POS)},
        };
    }

    @Test(dataProvider = "multiTileData")
    public void multiTileDataTest(final List<Integer> tiles, final Integer startingTileIndex, final List<File> files) {
        final IlluminaFileMap fm = new IlluminaFileMap();
        for(int i = 0; i < tiles.size(); i++) {
            fm.put(tiles.get(i), files.get(i));
        }

        final PosParser parser = (startingTileIndex == null) ? new PosParser(fm, IlluminaFileUtil.SupportedIlluminaFormat.Pos) :
                                                               new PosParser(fm, tiles.get(startingTileIndex), IlluminaFileUtil.SupportedIlluminaFormat.Pos);
        final List<TestResult> expectedResultsList = new ArrayList<TestResult>();
        final int t1 = (startingTileIndex != null) ? startingTileIndex : 0;
        for(int i = t1; i < tiles.size(); i++) {
            expectedResultsList.addAll(TEST_DATA.get(files.get(i).getName()));
        }

        final Iterator<TestResult> expected = expectedResultsList.iterator();

        int index = 0;
        while(expected.hasNext()) {
            Assert.assertTrue(parser.hasNext());
            compareResults(expected.next(), parser.next(), index);
            ++index;
        }

        Assert.assertFalse(parser.hasNext());
        parser.close();
    }
}
