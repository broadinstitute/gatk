package org.broadinstitute.hellbender.tools.picard.illumina.parser;

import htsjdk.samtools.util.CloseableIterator;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.*;

import static htsjdk.samtools.util.CollectionUtil.makeList;

public class PerTileParserTest {

    private static final Map<String, List<Integer>> FILE_TO_VALUE = new HashMap<>();
    private static final IlluminaFileMap FILE_MAP = new IlluminaFileMap();
    static {
        FILE_MAP.put(1, new File("s_1_1"));
        FILE_MAP.put(2, new File("s_1_2"));
        FILE_MAP.put(3, new File("s_1_3"));
        FILE_MAP.put(4, new File("s_1_4"));
        FILE_MAP.put(5, new File("s_1_5"));
        FILE_TO_VALUE.put("s_1_1", makeList(1, 2, 3, 4, 5));
        FILE_TO_VALUE.put("s_1_2", makeList(6, 7, 8, 9, 10));
        FILE_TO_VALUE.put("s_1_3", makeList(11, 12, 13, 14, 15));
        FILE_TO_VALUE.put("s_1_4", makeList(16, 17, 18, 19, 20));
        FILE_TO_VALUE.put("s_1_5", makeList(21, 22, 23, 24, 25));
    }

    class DummyDt implements IlluminaData {
        public DummyDt(final Integer value) {
            this.value = value;
        }
        public final Integer value;
    }

    class MockPerTileParser extends PerTileParser<DummyDt> {

        public MockPerTileParser(final IlluminaFileMap tilesToFiles) {
            super(tilesToFiles);
        }

        @Override
        protected CloseableIterator<DummyDt> makeTileIterator(final File file) {
            return new CloseableIterator<DummyDt>() {
                private final Iterator<Integer> values = FILE_TO_VALUE.get(file.getName()).iterator();

                @Override
                public void close() {

                }

                @Override
                public boolean hasNext() {
                    return values.hasNext();
                }

                @Override
                public DummyDt next() {
                    return new DummyDt(values.next());
                }

                @Override
                public void remove() {
                    throw new UnsupportedOperationException();
                }
            };
        }

        @Override
        public Set<IlluminaDataType> supportedTypes() {
            return null;  //To change body of implemented methods use File | Settings | File Templates.
        }
    }

    @Test
    public void basicIterationTest() {
        final IlluminaFileMap fm = new IlluminaFileMap();
        fm.put(1, new File("s_1_1"));
        fm.put(2, new File("s_1_2"));
        fm.put(3, new File("s_1_3"));
        fm.put(4, new File("s_1_4"));
        fm.put(5, new File("s_1_5"));

        final PerTileParser<DummyDt> ddts = new MockPerTileParser(fm);

        for(int i = 0; i < 25; i++) {
            Assert.assertTrue(ddts.hasNext());
            Assert.assertEquals(ddts.next().value, new Integer(i + 1));
        }

        Assert.assertFalse(ddts.hasNext());
    }

    @DataProvider(name="seekingTests")
    public Object[][] seekingTests() {
        return new Object[][] {
            {1,  4, null, null},
            {15, 1, null, null},
            {25, 3, null, null},
            {24, 5, null, null},
            {1,  3,   10, 1},
            {1,  3,   15, 2},
            {12, 2,   15, 4},
            {6,  3,   12, 5},
            {14, 5,   25, 2}
        };
    }

    @Test(dataProvider = "seekingTests")
    public void seekToTileTest(Integer firstSeekPos, Integer firstTile, Integer secondSeekPos, Integer secondTile) {
        final PerTileParser<DummyDt> ddts = new MockPerTileParser(FILE_MAP);

        for(int i = 1; i <= firstSeekPos; i++) {
            Assert.assertTrue(ddts.hasNext());
            Assert.assertEquals(ddts.next().value, new Integer(i));
        }

        ddts.seekToTile(firstTile);

        int start = firstTile * 5 - 4;
        if(secondSeekPos != null) {
            for(int i = start; i <= secondSeekPos; i++) {
                Assert.assertTrue(ddts.hasNext());
                Assert.assertEquals(ddts.next().value, new Integer(i));
            }
            ddts.seekToTile(secondTile);
            start = secondTile * 5 - 4;
        }

        for(int i = start; i <= 25; i++) {
            Assert.assertTrue(ddts.hasNext());
            Assert.assertEquals(ddts.next().value, new Integer(i));
        }

        Assert.assertFalse(ddts.hasNext());
    }

    @DataProvider(name="missingTiles")
    public Object[][] missingTiles() {
        return new Object[][] {
           {-1}, {10}, {Integer.MAX_VALUE}, {Integer.MIN_VALUE}
        };
    }

    @Test(expectedExceptions = IlluminaParserException.class, dataProvider="missingTiles")
    public void missingTileTest(final Integer missingTile) {
        final PerTileParser<DummyDt> ddts = new MockPerTileParser(FILE_MAP);
        ddts.seekToTile(missingTile);
    }

    @Test(expectedExceptions = IlluminaParserException.class)
    public void failVerifyTestTooManyTiles() {
        final PerTileParser<DummyDt> ddts = new MockPerTileParser(FILE_MAP);
        ddts.verifyData(makeList(1,2,3,4,5,6), null);
    }
    @Test(expectedExceptions = IlluminaParserException.class)
    public void failVerifyTestMissingTile() {
        final PerTileParser<DummyDt> ddts = new MockPerTileParser(FILE_MAP);
        ddts.verifyData(makeList(1,2,4,5), null);
    }

    @Test
    public void passVerifyTest() {
        final PerTileParser<DummyDt> ddts = new MockPerTileParser(FILE_MAP);
        ddts.verifyData(makeList(1,2,3,4,5), null);
    }
}
