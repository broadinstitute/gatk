package org.broadinstitute.hellbender.tools.picard.illumina.parser;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.*;

import static htsjdk.samtools.util.CollectionUtil.makeList;

public final class PerTilePerCycleParserTest {
    public static final List<Integer> DEFAULT_TILES = makeList(1, 2, 3, 4);
    public static final int[] DEFAULT_OUTPUT_LENGTHS = new int[]{10, 5, 5};
    public static final int[] CYCLES = new int[]{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20};
    public static final int MAX_CYCLE = 20;
    public static final OutputMapping OUTPUT_MAPPING = new OutputMapping(new ReadStructure("20T"));

    private class MockCycledIlluminaData implements IlluminaData {
        private final List<String> values;

        public MockCycledIlluminaData() {
            this.values = new ArrayList<String>();
        }

        public void addValue(final String value) {
            values.add(value);
        }

        public List<String> getValues() {
            return values;
        }
    }

    class MockPerTileCycleParser extends PerTileCycleParser<MockCycledIlluminaData> {
        private final int[] expectedOutputLengths;

        public MockPerTileCycleParser(final File directory, final int lane, final CycleIlluminaFileMap tilesToCycleFiles, final OutputMapping outputMapping) {
            super(directory, lane, tilesToCycleFiles, outputMapping);
            expectedOutputLengths = outputMapping.getOutputReadLengths();
            this.initialize();
        }

        @Override
        protected CycleFilesParser<MockCycledIlluminaData> makeCycleFileParser(final List<File> files) {
            return new CycleFilesParser<MockCycledIlluminaData>() {
                int currentCycle = 0;

                @Override
                public void close() {

                }

                @Override
                public MockCycledIlluminaData next() {
                    final MockCycledIlluminaData ild = new MockCycledIlluminaData();
                    if (!hasNext()) {
                        throw new NoSuchElementException();
                    }

                    ild.addValue(str_del(files.get(currentCycle++)));
                    return ild;
                }

                @Override
                public boolean hasNext() {
                    return currentCycle < MAX_CYCLE;
                }
            };
        }

        @Override
        public void initialize() {
            seekToTile(currentTile);
        }

        public Set<IlluminaDataType> supportedTypes() {
            return null;
        }

        @Override
        public void close() {
            //no-op
        }
    }

    public List<String> getFileNames(final List<Integer> tiles) {
        final List<String> fileNames = new ArrayList<String>();
        for (final Integer tile : tiles) {
            for (int i = 1; i <= MAX_CYCLE; i++) {
                fileNames.add(str_del(tile, i));
            }
        }
        return fileNames;
    }

    public CycleIlluminaFileMap getIlluminaFileMaps(final List<Integer> tiles, final int[] cycles) {
        final CycleIlluminaFileMap cycleFileMap = new CycleIlluminaFileMap();
        for (final int cycle : cycles) {
            final IlluminaFileMap fileMap = new IlluminaFileMap();
            for (final Integer tile : tiles) {
                fileMap.put(tile, new File(str_del(tile, cycle)));
            }
            cycleFileMap.put(cycle, fileMap);
        }
        return cycleFileMap;
    }

    public static String str_del(final Object... objects) {
        String out = objects[0].toString();
        for (int i = 1; i < objects.length; i++) {
            out += "_" + objects[i];
        }
        return out;
    }

    @Test
    public void basicIterationTest() {
        final List<String> expectedValues = getFileNames(DEFAULT_TILES);
        final PerTileCycleParser<MockCycledIlluminaData> parser = makeParser();

        int index = 0;
        while (parser.hasNext()) {
            index = compareValues(parser.next().values, expectedValues, index);
        }

        Assert.assertEquals(index, expectedValues.size());
    }


    private int compareValues(final List<String> parserValues, final List<String> expectedValues, int index) {
        for (final String parserValue : parserValues) {
            Assert.assertTrue(index < expectedValues.size());
            Assert.assertEquals(parserValue, expectedValues.get(index), "With index " + index);
            ++index;
        }

        return index;
    }

    public PerTileCycleParser<MockCycledIlluminaData> makeParser() {
        final CycleIlluminaFileMap fileMap = getIlluminaFileMaps(DEFAULT_TILES, CYCLES);
        return new MockPerTileCycleParser(new File("FakeFile"), 1, fileMap, OUTPUT_MAPPING);
    }

    @DataProvider(name = "seekingTests")
    public Object[][] seekingTests() {
        return new Object[][]{
                {1, 3, null, null},
                {22, 1, null, null},
                {38, 2, null, null},
                {75, 4, null, null},
                {1, 3, 70, 1},
                {1, 3, 45, 2},
                {12, 2, 59, 4},
                {45, 3, 70, 3},
                {14, 1, 5, 2}
        };
    }


    @Test(dataProvider = "seekingTests")
    public void seekingIterationTest(final Integer seekPos1, final Integer newTile1, final Integer seekPos2, final Integer newTile2) {
        final List<String> expectedValues = getFileNames(DEFAULT_TILES);
        final PerTileCycleParser<MockCycledIlluminaData> parser = makeParser();

        int index = 0;
        for (int i = 0; i <= seekPos1; i++) {
            Assert.assertTrue(parser.hasNext());
            index = compareValues(parser.next().values, expectedValues, index);
        }

        parser.seekToTile(newTile1);

        index = (newTile1 - 1) * MAX_CYCLE;
        if (seekPos2 != null) {
            for (int i = index; i <= seekPos2; i++) {
                Assert.assertTrue(parser.hasNext());
                index = compareValues(parser.next().values, expectedValues, index);
            }

            parser.seekToTile(newTile2);
            index = (newTile2 - 1) * MAX_CYCLE;
        }

        for (int i = index; i < MAX_CYCLE * DEFAULT_TILES.size(); i++) {
            Assert.assertTrue(parser.hasNext());
            index = compareValues(parser.next().values, expectedValues, index);
        }

        Assert.assertFalse(parser.hasNext());

    }
}
