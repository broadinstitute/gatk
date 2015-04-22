package org.broadinstitute.hellbender.tools.picard.illumina.parser;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.lang.reflect.Array;
import java.util.*;

public final class FilterParserTest {
    private static final File TestDataDir = new File("src/test/resources/org/broadinstitute/hellbender/tools/picard/illumina/parserTests");
    private static final File BaseCallsDir = new File(TestDataDir, "filterParser");

    private static final boolean F = false;
    private static final boolean T = true;
    public static final Boolean[] s_1_0001_filter = {F, F, T, F, T, T, F, F, T, F, T, T, F, T, T};
    public static final Boolean[] s_1_0002_filter = {F, F, F, F, F, T, T, F, F, F, F, T, F, F, F, T, T, T, T, T};
    public static final Boolean[] s_1_0003_filter = {T, F, F, F, T, T, F, F};
    public static final Boolean[] s_1_0004_filter = {T, F, F, F, F, T, F, F, F, F, T, F, F, F, F, T, T, F, F, T,
            T, F, T, F, T, T, F, F, T, T, F, T, T, T, T, F, T, T, T, T,
            F, T, F, F, F, T, F, T, F, T, T, F, F, T, T, T, F, F, T, T};

    public static final Boolean[] tileToValue(int tile) {
        switch (tile) {
            case 1:
                return s_1_0001_filter;
            case 2:
                return s_1_0002_filter;
            case 3:
                return s_1_0003_filter;
            case 4:
                return s_1_0004_filter;
        }

        throw new RuntimeException("You shouldn't reach this statement!");
    }

    public static final Boolean[] allTilesToValues(Integer[] tiles) {
        Boolean[][] values = new Boolean[tiles.length][];
        for (int i = 0; i < tiles.length; i++) {
            values[i] = tileToValue(tiles[i]);
        }

        return arrayFlatten(values);
    }

    public static final List<Integer> arrayToList(final Integer[] array) {
        final List<Integer> list = new ArrayList<Integer>();
        for (int item : array) {
            list.add(item);
        }
        return list;
    }

    public static final <T> T[] arrayFlatten(final T[][] arrays) {
        int total = 0;
        for (T[] arr : arrays) {
            total += arr.length;
        }

        int resultIndex = 0;
        final T[] result = (T[]) Array.newInstance(arrays[0][0].getClass(), total);
        for (int i = 0; i < arrays.length; i++) {
            System.arraycopy(arrays[i], 0, result, resultIndex, arrays[i].length);
            resultIndex += arrays[i].length;
        }
        return result;
    }

    @DataProvider(name = "passingTiles")
    public Object[][] passingTiles() {
        return new Object[][]{
                {new Integer[]{1}},
                {new Integer[]{2}},
                {new Integer[]{4}},
                {new Integer[]{1, 4}},
                {new Integer[]{2, 3}},
                {new Integer[]{1, 2, 3, 4}}
        };
    }

    @Test(dataProvider = "passingTiles")
    public void passingParserTest(Integer[] tiles) {
        final IlluminaFileUtil fUtil = new IlluminaFileUtil(BaseCallsDir, 1);

        final Boolean[] values = allTilesToValues(tiles);
        final List<Integer> tileList = arrayToList(tiles);

        final FilterParser fp = new FilterParser(((PerTileFileUtil) fUtil.getUtil(IlluminaFileUtil.SupportedIlluminaFormat.Filter)).getFiles(tileList));
        fp.verifyData(tileList, null); //filter parser doesn't care about cycle

        for (final boolean expectedValue : values) {
            Assert.assertTrue(fp.hasNext());
            Assert.assertEquals(fp.next().isPf(), expectedValue);
        }

        Assert.assertFalse(fp.hasNext());

        //seek back to the beginning and do it again!
        fp.seekToTile(tiles[0]);
        for (final boolean expectedValue : values) {
            Assert.assertTrue(fp.hasNext());
            Assert.assertEquals(fp.next().isPf(), expectedValue);
        }

        Assert.assertFalse(fp.hasNext());
        fp.close();
    }

    @DataProvider(name = "seekToTile")
    public Object[][] seekToTile() {
        return new Object[][]{
                //seek after how many, tiles to load, expected values
                {0, 2, new Integer[]{2}, new Boolean[][]{s_1_0002_filter}},
                {4, 4, new Integer[]{1, 4}, new Boolean[][]{new Boolean[]{F, F, T, F}, s_1_0004_filter}},
                {0, 3, new Integer[]{2, 3}, new Boolean[][]{s_1_0003_filter}},
                {15, 3, new Integer[]{1, 2, 3, 4}, new Boolean[][]{s_1_0001_filter, s_1_0003_filter, s_1_0004_filter}}
        };
    }

    @Test(dataProvider = "seekToTile")
    public void passingSeekingParserTest(int skipBefore, int tileToSkipTo, Integer[] tiles, Boolean[][] expectedValues) {
        final IlluminaFileUtil fUtil = new IlluminaFileUtil(BaseCallsDir, 1);

        final Boolean[] values = arrayFlatten(expectedValues);
        final List<Integer> tileList = arrayToList(tiles);

        final FilterParser fp = new FilterParser(((PerTileFileUtil) fUtil.getUtil(IlluminaFileUtil.SupportedIlluminaFormat.Filter)).getFiles(tileList));

        int read = 0;
        for (final boolean expectedValue : values) {
            if (read == skipBefore) {
                fp.seekToTile(tileToSkipTo);
            }

            Assert.assertTrue(fp.hasNext());
            Assert.assertEquals(fp.next().isPf(), expectedValue, " Failed on read: " + read);
            ++read;
        }

        Assert.assertFalse(fp.hasNext());
        fp.close();
    }

    @DataProvider(name = "failingVerifyTiles")
    public Object[][] failingVerifyTiles() {
        return new Object[][]{
                {new Integer[]{1}, new Integer[]{2}},
                {new Integer[]{2}, new Integer[]{1}},
                {new Integer[]{4}, new Integer[]{5}},
                {new Integer[]{1, 4}, new Integer[]{1, 3, 4}},
                {new Integer[]{2, 3}, new Integer[]{2, 3, 4}},
                {new Integer[]{1, 2, 3, 4}, new Integer[]{1, 2, 3, 4, 5}},
                {new Integer[]{2, 3, 4}, new Integer[]{1, 2, 3, 4}}
        };
    }

    @Test(dataProvider = "failingVerifyTiles", expectedExceptions = IlluminaParserException.class)
    public void verifyDataTest(final Integer[] initTiles, final Integer[] verifyTiles) {
        final IlluminaFileUtil fUtil = new IlluminaFileUtil(BaseCallsDir, 1);
        final List<Integer> initTileList = arrayToList(initTiles);
        final List<Integer> verifyTileList = arrayToList(verifyTiles);

        final FilterParser fp = new FilterParser(((PerTileFileUtil) fUtil.getUtil(IlluminaFileUtil.SupportedIlluminaFormat.Filter)).getFiles(initTileList));
        fp.verifyData(verifyTileList, null); //filter parser doesn't care about cycle values
        fp.close();
    }
}
