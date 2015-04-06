package org.broadinstitute.hellbender.tools.picard.illumina.parser;

import htsjdk.samtools.util.IOUtil;
import org.testng.Assert;
import org.testng.annotations.AfterMethod;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import org.broadinstitute.hellbender.tools.picard.illumina.parser.IlluminaFileUtil.SupportedIlluminaFormat;

import java.io.DataOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import static htsjdk.samtools.util.CollectionUtil.makeList;

public class IlluminaFileUtilTest {
    private static final int DEFAULT_LANE = 7;
    private static final List<Integer> DEFAULT_TILES = makeList(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12);
    private static final List<Integer> DEFAULT_TILE_TEST_SUBSET = makeList(1, 4, 5, 6, 9, 10);
    private static final int[] DEFAULT_CYCLES = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20};
    private static final int DEFAULT_LAST_CYCLE = 20;

    // We have data for these that should agree
    private static final List<SupportedIlluminaFormat> FORMATS_TO_TEST = Arrays.asList(
            SupportedIlluminaFormat.Bcl,
            SupportedIlluminaFormat.Locs,
            SupportedIlluminaFormat.Clocs,
            SupportedIlluminaFormat.Pos,
            SupportedIlluminaFormat.Filter,
            SupportedIlluminaFormat.Barcode);


    private File intensityDir;
    private File basecallDir;

    @BeforeMethod
    private void setUp() throws Exception {
        intensityDir = IOUtil.createTempDir("ift_test", "Intensities");
        basecallDir = new File(intensityDir, "BaseCalls");
        if (!basecallDir.mkdir()) {
            throw new RuntimeException("Couldn't make basecalls dir " + basecallDir.getAbsolutePath());
        }
    }

    @AfterMethod
    private void tearDown() {
        IOUtil.deleteDirectoryTree(intensityDir);
    }

    @DataProvider(name = "validLanes")
    public Object[][] validLanes() {
        return new Object[][]{
                {0, "s_0_1111.test"},
                {1, "s_1_23.test"},
                {10, "s_10_1.test"}
        };
    }

    public void regexMatches(final String regex, final String toMatch) {
        regexMatches(regex, toMatch, true);
    }

    public void regexMatches(final String regex, final String toMatch, final boolean expectedResult) {
        final Pattern pt = Pattern.compile(regex);
        final Matcher ma = pt.matcher(toMatch);
        Assert.assertEquals(ma.matches(), expectedResult);
    }

    @Test(dataProvider = "validLanes")
    public void regexTests(final int lane, final String ltExample) {
        regexMatches(ParameterizedFileUtil.makeLaneTileRegex(".test", lane), ltExample);
    }

    @DataProvider(name = "validLanesInvalidRegexes")
    public Object[][] validLanesInvalidRegexes() {
        return new Object[][]{
                {0, "s_-0_111"},
                {1, "s_1_A3"},
                {10, "s_-100_1"},
                {20, "s_21_1"}
        };
    }

    @Test(dataProvider = "validLanesInvalidRegexes")
    public void notMatchingRegexTest(final int lane, final String ltExample) {
        regexMatches(ParameterizedFileUtil.makeLaneTileRegex(".test", lane) , ltExample, false);
    }

    @DataProvider(name = "invalidLanes")
    public Object[][] invalidLanes() {
        return new Object[][]{
                {-1000},
                {-10},
                {-1}
        };
    }

    @Test(dataProvider = "invalidLanes", expectedExceptions = IlluminaParserException.class)
    public void invalidLaneForLTRegex(final int lane) {
        ParameterizedFileUtil.makeLaneTileRegex(".test", lane);
    }

    public void assertDefaults(final IlluminaFileUtil fileUtil, final Integer lane, final List<SupportedIlluminaFormat> formatsToTest) {
        if (lane == null) {
            Assert.assertEquals(fileUtil.getLane(), DEFAULT_LANE);
        } else {
            Assert.assertEquals(new Integer(fileUtil.getLane()), lane);
        }

        Assert.assertEquals(fileUtil.getUtil(SupportedIlluminaFormat.Barcode).getTiles(), DEFAULT_TILES);

        Assert.assertEquals(fileUtil.getUtil(SupportedIlluminaFormat.Bcl).getTiles(), DEFAULT_TILES);

        Assert.assertEquals(fileUtil.getUtil(SupportedIlluminaFormat.Pos).getTiles(), DEFAULT_TILES);

        Assert.assertEquals(fileUtil.getUtil(SupportedIlluminaFormat.Locs).getTiles(), DEFAULT_TILES);

        Assert.assertEquals(fileUtil.getUtil(SupportedIlluminaFormat.Clocs).getTiles(), DEFAULT_TILES);

        Assert.assertEquals(fileUtil.getUtil(SupportedIlluminaFormat.Filter).getTiles(), DEFAULT_TILES);

        final Set<Integer> detectedCycles = ((PerTilePerCycleFileUtil) fileUtil.getUtil(SupportedIlluminaFormat.Bcl)).getDetectedCycles();
        Assert.assertEquals(detectedCycles.size(), DEFAULT_CYCLES.length);
        int i = 0;
        for(final Integer cycle : detectedCycles){
            Assert.assertEquals(cycle.intValue(), DEFAULT_CYCLES[i++], "Elements differ at index " + i);
        }

        Assert.assertEquals(fileUtil.getActualTiles(formatsToTest), DEFAULT_TILES);
    }

    @Test
    public void passNewUtilTest() {
        for (final SupportedIlluminaFormat format : SupportedIlluminaFormat.values()) {
            makeFiles(format, intensityDir, DEFAULT_LANE, DEFAULT_TILES, DEFAULT_CYCLES);
            makeFiles(format, intensityDir, DEFAULT_LANE + 1, DEFAULT_TILES, DEFAULT_CYCLES, ".gz");
            makeFiles(format, intensityDir, DEFAULT_LANE + 2, DEFAULT_TILES, DEFAULT_CYCLES, ".bz2");
        }

        final Set<SupportedIlluminaFormat> formatsToTest = new HashSet<SupportedIlluminaFormat>();
        // TODO: I can't be bothered to build files for these.  AW
        Collections.addAll(formatsToTest, SupportedIlluminaFormat.values());
        formatsToTest.remove(SupportedIlluminaFormat.MultiTileBcl);
        formatsToTest.remove(SupportedIlluminaFormat.MultiTileFilter);
        formatsToTest.remove(SupportedIlluminaFormat.MultiTileLocs);
        final ArrayList<SupportedIlluminaFormat> formatsList = new ArrayList<SupportedIlluminaFormat>(formatsToTest);

        for (int i = 0; i < 3; i++) {
            final IlluminaFileUtil fileUtil = new IlluminaFileUtil(new File(intensityDir, "BaseCalls"), DEFAULT_LANE + i);
            Assert.assertEquals(fileUtil.getActualTiles(formatsList), DEFAULT_TILES);
            assertDefaults(fileUtil, DEFAULT_LANE + i, formatsList);
        }
    }

    @Test
    public void passingVerifyTest() {
        for (final SupportedIlluminaFormat format : SupportedIlluminaFormat.values()) {
            makeFiles(format, intensityDir, DEFAULT_LANE, DEFAULT_TILES, DEFAULT_CYCLES);
            makeFiles(format, intensityDir, DEFAULT_LANE + 1, DEFAULT_TILES, DEFAULT_CYCLES, ".gz");
            makeFiles(format, intensityDir, DEFAULT_LANE + 2, DEFAULT_TILES, DEFAULT_CYCLES, ".bz2");
        }

        for (int i = 0; i < 3; i++) {
            final IlluminaFileUtil fileUtil = new IlluminaFileUtil(new File(intensityDir, "BaseCalls"), DEFAULT_LANE + i);


            for (final SupportedIlluminaFormat format : FORMATS_TO_TEST) {
                Assert.assertEquals(new ArrayList<String>(), fileUtil.getUtil(format).verify(DEFAULT_TILES, DEFAULT_CYCLES));
            }
        }
    }

    @DataProvider(name = "missingTileFormats")
    public Object[][] missingTileFormats() {
        return new Object[][]{
                {
                        1,
                        makeList(SupportedIlluminaFormat.Bcl, SupportedIlluminaFormat.Barcode),
                        makeList(SupportedIlluminaFormat.Bcl, SupportedIlluminaFormat.Barcode),
                        makeList("BaseCalls/s_1_0007_barcode.txt.gz"),
                        ".gz"
                },

                {
                        2,
                        Arrays.asList(SupportedIlluminaFormat.values()),
                        Arrays.asList(SupportedIlluminaFormat.values()),
                        makeCycleFileList(new File("BaseCalls"), ".bcl", 2, DEFAULT_CYCLES, 2),
                        ".gz"
                },
                {
                        3,
                        Arrays.asList(SupportedIlluminaFormat.values()),
                        Arrays.asList(SupportedIlluminaFormat.values()),
                        makeList("BaseCalls/L003/C1.1/s_3_2.bcl"),
                        ".bz2"
                },
                {
                        4,
                        Arrays.asList(SupportedIlluminaFormat.values()),
                        Arrays.asList(SupportedIlluminaFormat.Pos, SupportedIlluminaFormat.Locs),
                        makeList("s_4_10_pos.txt", "L004/s_4_2.locs"),
                        null
                }
        };
    }

    public static void emptyRelativeFiles(final File baseFile, final List<String> relativeFilesToDelete) {
        for (final String relativeFile : relativeFilesToDelete) {
            final File actualFile = new File(baseFile, relativeFile);


            if (!actualFile.exists()) {
                throw new RuntimeException("Trying to empty a non-existent file" + actualFile.getAbsolutePath());
            }

            if (actualFile.isDirectory()) {
                throw new RuntimeException("Trying to empty a directory(" + actualFile.getAbsolutePath() + ")");
            } else {
                if (!actualFile.delete()) {
                    throw new RuntimeException("Couldn't remove previous file when emptying(" + actualFile.getAbsolutePath() + ")");
                } else {
                    try {
                        if (!actualFile.createNewFile()) {
                            throw new RuntimeException("Couldn't create empty file: " + actualFile.getAbsolutePath() + ")");
                        }
                    } catch (final IOException ioe) {
                        throw new RuntimeException(ioe);
                    }
                }
            }
            if (!actualFile.exists()) {
                throw new RuntimeException("File should exist: " + actualFile);
            }
        }
    }

    public static void deleteRelativeFiles(final File baseFile, final List<String> relativeFilesToDelete) {
        for (final String relativeFile : relativeFilesToDelete) {
            final File actualFile = new File(baseFile, relativeFile);


            if (!actualFile.exists()) {
                throw new RuntimeException("Trying to delete a non-existent file" + actualFile.getAbsolutePath());
            }

            if (actualFile.isDirectory()) {
                IOUtil.deleteDirectoryTree(actualFile);
            } else {
                IOUtil.deleteFiles(actualFile);
            }
            if (actualFile.exists()) {
                throw new RuntimeException("File still exists after calling delete: " + actualFile);
            }
        }
    }

    public final void deleteRelativeFiles(final List<String> relativeFilesToDelete) {
        deleteRelativeFiles(intensityDir, relativeFilesToDelete);
    }

    @Test(dataProvider = "missingTileFormats")
    public void missingTileTest(final int lane,
                                final List<SupportedIlluminaFormat> formats,
                                final List<SupportedIlluminaFormat> formatsToGetTiles,
                                final List<String> relativeFilesToDelete,
                                final String compression) {
        for (final SupportedIlluminaFormat format : formats) {
            makeFiles(format, intensityDir, lane, DEFAULT_TILES, DEFAULT_CYCLES, compression);
        }

        deleteRelativeFiles(relativeFilesToDelete);

        IlluminaParserException ipe = null;
        try {
            final IlluminaFileUtil fUtil = new IlluminaFileUtil(new File(intensityDir, "BaseCalls"), lane);
            fUtil.getActualTiles(formatsToGetTiles);
        } catch (final IlluminaParserException exception) {
            ipe = exception;
        }

        Assert.assertNotNull(ipe, "Didn't raise a Picard Exception for missing tile!");
        Assert.assertTrue(ipe.getMessage().contains("Formats do not have the same number of tiles! "), "Wrong exception thrown for missing tile!");
    }

    @DataProvider(name = "perTileFileFormats")
    public Object[][] perTileFileUtils() {
        return new Object[][]{
                {SupportedIlluminaFormat.Locs, null, false, laneDir(DEFAULT_LANE)},
                {SupportedIlluminaFormat.Clocs, null, false, laneDir(DEFAULT_LANE)},
                {SupportedIlluminaFormat.Pos, ".gz", false, null},
                {SupportedIlluminaFormat.Pos, null, false, null},
                {SupportedIlluminaFormat.Filter, null, true, "BaseCalls/" + laneDir(DEFAULT_LANE)},
                {SupportedIlluminaFormat.Barcode, ".bz2", true, "BaseCalls"}
        };
    }

    public File makePerTileFile(final File parentDir, final int lane, final int tile, final String extension, final String compression, final boolean longFormat) {
        return new File(parentDir, "s_" + lane + "_" + longTile(tile, longFormat) + extension + (compression != null ? compression : ""));
    }

    public void testDefaultPerTileUtil(final PerTileFileUtil ptfu, final String compression, final boolean longFormat, final File parentDir) {
        final IlluminaFileMap fm = ptfu.getFiles();
        final IlluminaFileMap fmWTiles = ptfu.getFiles(DEFAULT_TILES);

        Assert.assertEquals(fm.size(), DEFAULT_TILES.size());

        for (final Integer tile : DEFAULT_TILES) {
            final File tFile = fm.get(tile);
            final File tFile2 = fmWTiles.get(tile);
            Assert.assertEquals(tFile.getAbsolutePath(), tFile2.getAbsolutePath());
            Assert.assertEquals(tFile, makePerTileFile(parentDir, DEFAULT_LANE, tile, ptfu.extension, compression, longFormat));
            Assert.assertTrue(tFile.exists());
            Assert.assertTrue(tFile.length() > 0);
        }

        final List<Integer> tiles = new ArrayList<Integer>(DEFAULT_TILE_TEST_SUBSET);
        final IlluminaFileMap subsetMap = ptfu.getFiles(DEFAULT_TILE_TEST_SUBSET);
        for (final Integer tile : subsetMap.keySet()) {
            tiles.remove(tile);
            Assert.assertTrue(DEFAULT_TILE_TEST_SUBSET.contains(tile));
            final File tFile = subsetMap.get(tile);
            Assert.assertEquals(tFile, makePerTileFile(parentDir, DEFAULT_LANE, tile, ptfu.extension, compression, longFormat));
            Assert.assertTrue(tFile.exists());
            Assert.assertTrue(tFile.length() > 0);
        }

        Assert.assertTrue(tiles.isEmpty());
    }

    @Test(dataProvider = "perTileFileFormats")
    public void perTileFileUtilsTest(final SupportedIlluminaFormat format, final String compression, final boolean longFormat, final String parentDir) {
        makeFiles(format, intensityDir, DEFAULT_LANE, DEFAULT_TILES, DEFAULT_CYCLES, compression);

        final IlluminaFileUtil fileUtil = new IlluminaFileUtil(basecallDir, DEFAULT_LANE);
        final PerTileFileUtil ptfu = (PerTileFileUtil) fileUtil.getUtil(format);

        Assert.assertTrue(ptfu.filesAvailable());
        testDefaultPerTileUtil(ptfu, compression, longFormat, (parentDir == null) ? intensityDir : new File(intensityDir, parentDir));

        final IlluminaFileUtil noFilesFu = new IlluminaFileUtil(basecallDir, DEFAULT_LANE + 20);
        final PerTileFileUtil noFilesPtfu = (PerTileFileUtil) noFilesFu.getUtil(format);
        Assert.assertFalse(noFilesPtfu.filesAvailable());
        Assert.assertTrue(noFilesPtfu.getFiles().isEmpty());
        Assert.assertTrue(noFilesPtfu.getFiles(DEFAULT_TILES).isEmpty());
    }

    public File makePerTilePerCycleFilePath(final File parentDir, final int lane, final int tile, final int cycle, final String extension) {
        return new File(parentDir, "C" + cycle + ".1/s_" + lane + "_" + tile + extension);
    }

    public void testDefaultPerTilePerCycleUtil(final PerTilePerCycleFileUtil pcfu, final File parentDir, final int[] cycles) {
        final CycleIlluminaFileMap cfm = pcfu.getFiles(cycles);
        final CycleIlluminaFileMap cfmWTiles = pcfu.getFiles(DEFAULT_TILES, cycles);
        final CycleIlluminaFileMap cfmNoCycles;
        if (Arrays.equals(cycles, DEFAULT_CYCLES)) {
            cfmNoCycles = pcfu.getFiles();
        } else {
            cfmNoCycles = null;
        }

        Assert.assertEquals(cfm.size(), cycles.length);

        for (final int cycle : cycles) {
            final IlluminaFileMap tFileIter = cfm.get(cycle);
            final IlluminaFileMap tFileIter2 = cfmWTiles.get(cycle);
            final IlluminaFileMap tFileIter3;
            if (cfmNoCycles != null) {
                tFileIter3 = cfmNoCycles.get(cycle);
            } else {
                tFileIter3 = null;
            }

            for (final Integer tile : DEFAULT_TILES) {
                final File tcFile = tFileIter.get(tile);
                final File tcFile2 = tFileIter2.get(tile);

                Assert.assertEquals(tcFile.getAbsolutePath(), tcFile2.getAbsolutePath());
                if (tFileIter3 != null) {
                    final File tfFile3 = tFileIter3.get(tile);
                    Assert.assertEquals(tcFile.getAbsolutePath(), tfFile3.getAbsolutePath());
                }

                Assert.assertEquals(tcFile, makePerTilePerCycleFilePath(parentDir, DEFAULT_LANE, tile, cycle, pcfu.extension));
                Assert.assertTrue(tcFile.exists());
                Assert.assertTrue(tcFile.length() > 0);
            }
        }
    }


    public void testSubsetDefaultPerTilePerCycleUtil(final PerTilePerCycleFileUtil pcfu, final File parentDir, final int[] cycles) {
        final List<Integer> tiles = new ArrayList<Integer>(DEFAULT_TILE_TEST_SUBSET);
        final CycleIlluminaFileMap subsetMap = pcfu.getFiles(DEFAULT_TILE_TEST_SUBSET, cycles);
        final CycleIlluminaFileMap cfmNoCycles;
        if (Arrays.equals(cycles, DEFAULT_CYCLES)) {
            cfmNoCycles = pcfu.getFiles(DEFAULT_TILE_TEST_SUBSET);
        } else {
            cfmNoCycles = null;
        }

        for (final int cycle : cycles) {
            final IlluminaFileMap tFileIter = subsetMap.get(cycle);
            final IlluminaFileMap tFileIter2;
            if (cfmNoCycles != null) {
                tFileIter2 = cfmNoCycles.get(cycle);
            } else {
                tFileIter2 = null;
            }


            for (final Integer tile : subsetMap.get(cycle).keySet()) {
                Assert.assertTrue(DEFAULT_TILE_TEST_SUBSET.contains(tile));
                tiles.remove(tile);
                final File tcFile = tFileIter.get(tile);
                if (tFileIter2 != null) {
                    Assert.assertEquals(tcFile, tFileIter2.get(tile));
                }
                Assert.assertEquals(tcFile, makePerTilePerCycleFilePath(parentDir, DEFAULT_LANE, tile, cycle, pcfu.extension));
                Assert.assertTrue(tcFile.exists());
                Assert.assertTrue(tcFile.length() > 0);
            }
        }

        Assert.assertTrue(tiles.isEmpty());
    }

    public static int[] cycleRange(final Range range) {
        return cycleRange(range.start, range.end);
    }

    public static int[] cycleRange(final int start, final int end) {
        final int[] cycles = new int[end - start + 1];
        for (int i = 0; i < cycles.length; i++) {
            cycles[i] = start + i;
        }

        return cycles;
    }

    public static int[] cycleRange(final int end) {
        return cycleRange(1, end);
    }

    @DataProvider(name = "perTilePerCycleFileFormats")
    public Object[][] perTilePerCycleFileFormats() {
        return new Object[][]{
                {SupportedIlluminaFormat.Bcl, "BaseCalls/" + laneDir(DEFAULT_LANE), DEFAULT_CYCLES, false, false},
                {SupportedIlluminaFormat.Bcl, "BaseCalls/" + laneDir(DEFAULT_LANE), cycleRange(4), true, true},
        };
    }

    @Test(dataProvider = "perTilePerCycleFileFormats")
    public void perTilePerCycleFileUtilsTest(final SupportedIlluminaFormat format, final String parentDir,
                                             final int[] cycles, final boolean createEarlySkippedCycles,
                                             final boolean createLateSkippedCycles) {
        if (createEarlySkippedCycles) {
            makeFiles(format, intensityDir, DEFAULT_LANE, DEFAULT_TILES, cycleRange(1, cycles[0]), null);
        }

        makeFiles(format, intensityDir, DEFAULT_LANE, DEFAULT_TILES, cycles, null);

        if (createLateSkippedCycles) {
            makeFiles(format, intensityDir, DEFAULT_LANE, DEFAULT_TILES, cycleRange(cycles[cycles.length - 1] + 1, DEFAULT_LAST_CYCLE), null);
        }

        final IlluminaFileUtil fileUtil = new IlluminaFileUtil(basecallDir, DEFAULT_LANE);
        final PerTilePerCycleFileUtil pcfu = (PerTilePerCycleFileUtil) fileUtil.getUtil(format);

        Assert.assertTrue(pcfu.filesAvailable());
        testDefaultPerTilePerCycleUtil(pcfu, (parentDir == null) ? intensityDir : new File(intensityDir, parentDir), cycles);
        testSubsetDefaultPerTilePerCycleUtil(pcfu, (parentDir == null) ? intensityDir : new File(intensityDir, parentDir), cycles);

        final IlluminaFileUtil noFilesFu = new IlluminaFileUtil(basecallDir, DEFAULT_LANE + 20);
        final PerTilePerCycleFileUtil noFilesPcfu = (PerTilePerCycleFileUtil) noFilesFu.getUtil(format);

        Assert.assertFalse(noFilesPcfu.filesAvailable());
        Assert.assertTrue(noFilesPcfu.getFiles().isEmpty());
        Assert.assertTrue(noFilesPcfu.getFiles(DEFAULT_TILES).isEmpty());
    }

    @DataProvider(name = "missingCycleDataRanges")
    public Object[][] missingCycleDataRanges() {
        return new Object[][]{
                {makeList(new Range(10, 15))},
                {makeList(new Range(9, 12), new Range(14, 15))}
        };
    }

    @Test(expectedExceptions = IlluminaParserException.class, dataProvider = "missingCycleDataRanges")
    public void perTilePerCycleFileUtilsMissingCycleTest(final List<Range> cycleRangesToMake) {
        final SupportedIlluminaFormat format = SupportedIlluminaFormat.Bcl;

        for (final Range range : cycleRangesToMake) {
            makeFiles(format, intensityDir, DEFAULT_LANE, DEFAULT_TILES, cycleRange(range), null);
        }

        final IlluminaFileUtil fileUtil = new IlluminaFileUtil(basecallDir, DEFAULT_LANE);
        final PerTilePerCycleFileUtil pcfu = (PerTilePerCycleFileUtil) fileUtil.getUtil(format);

        Assert.assertTrue(pcfu.filesAvailable());
        final int[] cycles = cycleRange(9, 16);
        final CycleIlluminaFileMap cfm = pcfu.getFiles(cycles);
        cfm.assertValid(DEFAULT_TILES, cycles);
    }


    public static void makeFiles(final SupportedIlluminaFormat format, final File intensityDir, final int lane,
                                 final List<Integer> tiles, final int[] cycles) {
        makeFiles(format, intensityDir, lane, tiles, cycles, null);
    }

    public static void makeFiles(final SupportedIlluminaFormat format, final File intensityDir, final int lane,
                                 final List<Integer> tiles, final int[] cycles, final String compression) {
        String laneDir = String.valueOf(lane);
        while (laneDir.length() < 3) {
            laneDir = "0" + laneDir;
        }
        laneDir = "L" + laneDir;


        final File basecallDir = new File(intensityDir, "BaseCalls");
        final File basecallLaneDir = new File(basecallDir, laneDir);
        final File intensityLaneDir = new File(intensityDir, laneDir);

        switch (format) {
            //per tile formats
            case Barcode:
                makePerTileFiles(basecallDir, lane, tiles, maybeAddExt("_barcode.txt", compression), true);
                break;

            case Pos:
                makePerTileFiles(intensityDir, lane, tiles, maybeAddExt("_pos.txt", compression), false);
                break;

            case Locs:
                makePerTileFiles(intensityLaneDir, lane, tiles, maybeAddExt(".locs", null), false);
                break;

            case Clocs:
                makePerTileFiles(intensityLaneDir, lane, tiles, maybeAddExt(".clocs", null), false);
                break;

            case Filter:
                makePerTileFiles(basecallLaneDir, lane, tiles, maybeAddExt(".filter", null), true);
                break;

            //per tile per cycle formats
            case Bcl:
                makePerTilePerCycleFiles(basecallLaneDir, lane, tiles, cycles, ".bcl");
                break;
        }
    }

    private static void makePerTileFiles(final File parentDir, final int lane, final List<Integer> tiles, final String ext, final boolean longName) {
        if (!parentDir.exists()) {
            if (!parentDir.mkdir()) {
                throw new RuntimeException("Couldn't create directory " + parentDir.getAbsolutePath());
            }
        }

        for (final Integer tile : tiles) {
            writeNonEmptyFile(new File(parentDir, "s_" + lane + "_" + longTile(tile, longName) + ext));
        }
    }

    private static void makePerTilePerCycleFiles(final File parentDir, final int lane, final List<Integer> tiles, final int[] cycles, final String ext) {
        if (!parentDir.exists()) {
            if (!parentDir.mkdir()) {
                throw new RuntimeException("Couldn't create directory " + parentDir.getAbsolutePath());
            }
        }

        for (final int cycle : cycles) {
            final File cycleDir = new File(parentDir, "C" + cycle + ".1");
            if (!cycleDir.exists()) {
                if (!cycleDir.mkdir()) {
                    throw new RuntimeException("Couldn't create directory " + cycleDir.getAbsolutePath());
                }
            }

            for (final Integer tile : tiles) {
                writeNonEmptyFile(new File(cycleDir, "s_" + lane + "_" + tile + ext));
            }
        }
    }

    private static List<String> makeCycleFileList(final File dir, final String ext, final int lane, final int[] cycles, final int... tiles) {
        return makeCycleFileList(dir, ext, lane, cycles, false, tiles);
    }

    private static List<String> makeCycleFileList(final File dir, final String ext, final int lane, final int[] cycles, final boolean longFmt, final int... tiles) {
        final List<String> files = new ArrayList<String>();
        final File laneDir = new File(dir, laneDir(lane));

        for (final int cycle : cycles) {
            final File cycleDir = new File(laneDir, "C" + cycle + ".1");
            for (final Integer tile : tiles) {
                files.add(cycleDir + "/s_" + lane + "_" + longTile(tile, longFmt) + ext);
            }
        }

        return files;
    }

    private static void writeNonEmptyFile(final File file) {
        try {
            final OutputStream outputStream = new DataOutputStream(new FileOutputStream(file));
            final int expectedLength = 10;
            outputStream.write(expectedLength);
            // The negative beginning index is to accommodate the header. Fancy. Ever so fancy.
            for (int i = -3; i < expectedLength; i++) outputStream.write(0x0);
            outputStream.close();
        } catch (final IOException e) {
            throw new RuntimeException("Exception trying to create non-empty file!", e);
        }
    }

    private static String laneDir(final int lane) {
        String ldir = String.valueOf(lane);
        while (ldir.length() < 3) {
            ldir = "0" + ldir;
        }
        return "L" + ldir;
    }

    private static String longTile(final int tile, final boolean makeLong) {
        if (makeLong) {
            String lt = String.valueOf(tile);
            while (lt.length() < 4) {
                lt = "0" + lt;
            }
            return lt;
        } else {
            return String.valueOf(tile);
        }
    }

    private static String maybeAddExt(final String fileExt, final String compressionExt) {
        if (compressionExt != null) {
            return fileExt + compressionExt;
        } else {
            return fileExt;
        }
    }
}
