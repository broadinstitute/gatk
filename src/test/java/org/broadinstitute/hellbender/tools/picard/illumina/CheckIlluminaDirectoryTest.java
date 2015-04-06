package org.broadinstitute.hellbender.tools.picard.illumina;

import htsjdk.samtools.SAMException;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.RuntimeIOException;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.tools.picard.illumina.parser.IlluminaDataType;
import org.broadinstitute.hellbender.tools.picard.illumina.parser.IlluminaFileUtil;
import org.broadinstitute.hellbender.tools.picard.illumina.parser.IlluminaFileUtilTest;
import org.testng.Assert;
import org.testng.annotations.AfterMethod;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.ByteOrder;
import java.nio.MappedByteBuffer;
import java.nio.channels.FileChannel;
import java.util.*;

import static org.broadinstitute.hellbender.tools.picard.illumina.parser.IlluminaDataType.BaseCalls;
import static org.broadinstitute.hellbender.tools.picard.illumina.parser.IlluminaDataType.Position;
import static org.broadinstitute.hellbender.tools.picard.illumina.parser.IlluminaFileUtil.SupportedIlluminaFormat;
import static org.broadinstitute.hellbender.tools.picard.illumina.parser.IlluminaFileUtil.SupportedIlluminaFormat.*;

import static htsjdk.samtools.util.CollectionUtil.makeList;

public class CheckIlluminaDirectoryTest extends CommandLineProgramTest {

    private File illuminaDir;
    private File dataDir;
    private File interopDir;
    private File intensityDir;
    private File basecallDir;

    public String getCommandLineProgramName() {
        return CheckIlluminaDirectory.class.getSimpleName();
    }

    @BeforeMethod
    private void setUp() throws Exception {
        illuminaDir = IOUtil.createTempDir("ift_test", "IlluminaDir");

        interopDir = new File(illuminaDir, "InterOp");
        if (!interopDir.exists() && !interopDir.mkdir()) {
            throw new RuntimeException("Couldn't make interop dir " + interopDir.getAbsolutePath());
        }

        dataDir = new File(illuminaDir, "Data");
        if (!dataDir.exists() && !dataDir.mkdir()) {
            throw new RuntimeException("Couldn't make data dir " + dataDir.getAbsolutePath());
        }

        intensityDir = new File(dataDir, "Intensities");
        if (!intensityDir.exists() && !intensityDir.mkdir()) {
            throw new RuntimeException("Couldn't make intensity dir " + intensityDir.getAbsolutePath());
        }

        basecallDir = new File(intensityDir, "BaseCalls");
        if (!basecallDir.exists() && !basecallDir.mkdir()) {
            throw new RuntimeException("Couldn't make basecalls dir " + basecallDir.getAbsolutePath());
        }
    }

    @AfterMethod
    private void tearDown() {
        IOUtil.deleteDirectoryTree(dataDir);
        IOUtil.deleteDirectoryTree(basecallDir);
        IOUtil.deleteDirectoryTree(intensityDir);
    }

    public void makeFiles(final SupportedIlluminaFormat[] formats, final int lane, final List<Integer> tiles,
                          final int[] cycles) {
        for (final IlluminaFileUtil.SupportedIlluminaFormat format : formats) {
            IlluminaFileUtilTest.makeFiles(format, intensityDir, lane, tiles, cycles);
        }
    }

    public List<String> makeCheckerArgs(final File basecallDir, final int lane, final String readStructure,
                                        final IlluminaDataType[] dataTypes, final List<Integer> filterTiles,
                                        final boolean makeFakeFiles, final boolean createSymLinks) {

        final List<String> dataTypeArgs = new ArrayList<>();
        dataTypeArgs.add("--B"); dataTypeArgs.add(basecallDir.getAbsolutePath());
        dataTypeArgs.add("--L"); dataTypeArgs.add(Integer.toString(lane));
        dataTypeArgs.add("--RS"); dataTypeArgs.add(readStructure);
        dataTypeArgs.add("--F"); dataTypeArgs.add(Boolean.toString(makeFakeFiles));
        dataTypeArgs.add("--X"); dataTypeArgs.add(Boolean.toString(createSymLinks));

        for (final IlluminaDataType dt : dataTypes) {
            dataTypeArgs.add("--DT"); dataTypeArgs.add(dt.toString());
        }

        for (final Integer ft : filterTiles) {
            dataTypeArgs.add("--T"); dataTypeArgs.add(Integer.toString(ft));
        }

        return dataTypeArgs;
    }

    public File writeTileMetricsOutFile(final Map<Integer, List<Integer>> lanesToTiles) {
        return writeTileMetricsOutFile(interopDir, (byte) 2, (byte) 10, lanesToTiles);
    }

    public File writeTileMetricsOutFile(final File interopDir, final byte versionNumber, final byte recordSize,
                                        final Map<Integer, List<Integer>> lanesToTiles) {
        final File tileMetricsOut = new File(interopDir, "TileMetricsOut.bin");
        if (!tileMetricsOut.exists()) {
            try {
                if (!tileMetricsOut.createNewFile()) {
                    throw new RuntimeIOException("Could not create tileMetricsOut file(" + tileMetricsOut.getAbsolutePath() + ")");
                }
            } catch (final IOException e) {
                throw new RuntimeIOException("IOException creating tileMetricsOut file (" + tileMetricsOut + ") for writing!", e);
            }
        }

        int totalEntries = 0;
        for (final Map.Entry<Integer, List<Integer>> l2t : lanesToTiles.entrySet()) {
            totalEntries += l2t.getValue().size();
        }

        final MappedByteBuffer buf;
        try {
            final RandomAccessFile raf = new RandomAccessFile(tileMetricsOut, "rw");
            final FileChannel channel = raf.getChannel();
            buf = channel.map(FileChannel.MapMode.READ_WRITE, 0, 2 + 10 * totalEntries);
            buf.order(ByteOrder.LITTLE_ENDIAN);

            buf.put(versionNumber);
            buf.put(recordSize);

            for (final int lane : lanesToTiles.keySet()) {
                for (final int tile : lanesToTiles.get(lane)) {
                    buf.putShort((short) lane);
                    buf.putShort((short) tile);
                    buf.putShort((short) 0);
                    buf.putFloat(0F);
                }
            }

            buf.force();
            CloserUtil.close(channel);
            CloserUtil.close(raf);
        } catch (final IOException e) {
            throw new RuntimeIOException("IOException writing tileMetricsOut file (" + tileMetricsOut + ")", e);
        }

        return tileMetricsOut;
    }

    public static Map<Integer, List<Integer>> makeMap(final List<Integer> lanes, final List<List<Integer>> tiles) {
        final Map<Integer, List<Integer>> map = new HashMap<Integer, List<Integer>>();

        if (lanes.size() != tiles.size()) {
            throw new IllegalArgumentException("Number of lanes (" + lanes + ") does not equal number of tiles!");
        }

        for (int i = 0; i < lanes.size(); i++) {
            map.put(lanes.get(i), tiles.get(i));
        }

        return map;
    }

    @DataProvider(name = "positiveTestData")
    public Object[][] positiveTestData() {
        return new Object[][]{
                {
                        new SupportedIlluminaFormat[]{Bcl, Locs, Pos, Filter},
                        new IlluminaDataType[]{BaseCalls, IlluminaDataType.QualityScores, IlluminaDataType.Position,
                                IlluminaDataType.PF},
                        3, makeList(1101, 1201, 1301, 2101, 2201, 2301), IlluminaFileUtilTest.cycleRange(1, 50),
                        "25T25T", new ArrayList<Integer>()
                },
                {
                        new SupportedIlluminaFormat[]{Bcl, Locs, Filter},
                        new IlluminaDataType[]{BaseCalls, IlluminaDataType.QualityScores, IlluminaDataType.Position,
                                IlluminaDataType.PF},
                        2, makeList(1101, 1201, 1301, 2101, 2201, 2301), IlluminaFileUtilTest.cycleRange(1, 50),
                        "8S15T8S", new ArrayList<Integer>()
                },
                {
                        new SupportedIlluminaFormat[]{Bcl, Filter},
                        new IlluminaDataType[]{BaseCalls, IlluminaDataType.QualityScores, IlluminaDataType.PF},
                        2, makeList(1101, 1201, 1301, 2101, 2201, 2301), IlluminaFileUtilTest.cycleRange(1, 152),
                        "68T8B68T", new ArrayList<Integer>()
                },
                {
                        new SupportedIlluminaFormat[]{Bcl, Pos, Filter},
                        new IlluminaDataType[]{BaseCalls, IlluminaDataType.QualityScores, IlluminaDataType.Position,
                                IlluminaDataType.PF},
                        5, makeList(1101, 1201, 1301, 2101, 2201, 2301), IlluminaFileUtilTest.cycleRange(1, 50),
                        "25T25T", new ArrayList<Integer>()
                },
                {
                        new SupportedIlluminaFormat[]{Bcl, Pos, Filter},
                        new IlluminaDataType[]{BaseCalls, IlluminaDataType.QualityScores, IlluminaDataType.Position,
                                IlluminaDataType.PF},
                        5, makeList(1101, 1201, 1301, 2101, 2201, 2301), IlluminaFileUtilTest.cycleRange(1, 50),
                        "25T25T", makeList(1301, 2101)
                }
        };
    }

    //Note: The positiveTest and negativeTests don't actually test Qseqs (the Qseq in the first test case above is there to make sure
    //BCLs are preferred over Qseqs)

    @Test(dataProvider = "positiveTestData")
    public void positiveTests(final IlluminaFileUtil.SupportedIlluminaFormat[] formats,
                              final IlluminaDataType[] dataTypes,
                              final int lane,
                              final List<Integer> tiles,
                              final int[] cycles,
                              final String readStructure,
                              final List<Integer> filterTiles) {
        makeFiles(formats, lane, tiles, cycles);
        writeTileMetricsOutFile(makeMap(makeList(lane - 1, lane + 1, lane),
                makeList(makeList(1, 2, 3), tiles, tiles)));

        final List<String> args = makeCheckerArgs(basecallDir, lane, readStructure, dataTypes, filterTiles, false, false);
        Assert.assertEquals(runCommandLine(args), 0); // assert that the number of found errors is zero
    }

    @DataProvider(name = "negativeTestData")
    public Object[][] negativeTestData() {
        return new Object[][]{
                { //Completely missing data types
                        new SupportedIlluminaFormat[]{Bcl, Filter},
                        new IlluminaDataType[]{BaseCalls, IlluminaDataType.QualityScores, IlluminaDataType.PF,
                                IlluminaDataType.Position, IlluminaDataType.Barcodes},
                        new ArrayList<String>(),
                        new ArrayList<String>(),
                        2, makeList(1101, 1201, 1301, 2101, 2201, 2301), IlluminaFileUtilTest.cycleRange(1, 152),
                        "68T8B68T",
                        2, new ArrayList<Integer>(), true
                },
                {
                        new SupportedIlluminaFormat[]{Bcl, Filter},
                        new IlluminaDataType[]{BaseCalls, IlluminaDataType.QualityScores, IlluminaDataType.PF},
                        makeList("BaseCalls/L002/C13.1/s_2_1201.bcl", "BaseCalls/L002/C13.1/s_2_2101.bcl"),
                        makeList("BaseCalls/L002/s_2_2101.filter"),
                        2, makeList(1101, 1201, 1301, 2101, 2201, 2301), IlluminaFileUtilTest.cycleRange(1, 20), "13T",
                        3, new ArrayList<Integer>(), true
                },
                {
                        new SupportedIlluminaFormat[]{Bcl, Filter},
                        new IlluminaDataType[]{BaseCalls, IlluminaDataType.QualityScores, IlluminaDataType.PF},
                        new ArrayList<String>(),
                        new ArrayList<String>(),
                        5, makeList(1101, 1201, 1301, 2101, 2201, 2301), IlluminaFileUtilTest.cycleRange(1, 152),
                        "250T",
                        98, new ArrayList<Integer>(), true
                },
                {
                        new SupportedIlluminaFormat[]{Bcl, Filter},
                        new IlluminaDataType[]{BaseCalls, IlluminaDataType.QualityScores, IlluminaDataType.PF},
                        new ArrayList<String>(),
                        new ArrayList<String>(),
                        5, makeList(1101, 1201, 1301, 2101, 2201, 2301), IlluminaFileUtilTest.cycleRange(1, 152),
                        "250T",
                        98, makeList(1301, 2201), true
                }
        };
    }

    @Test(dataProvider = "negativeTestData")
    public void negativeTests(final IlluminaFileUtil.SupportedIlluminaFormat[] formats,
                              final IlluminaDataType[] dataTypes,
                              final List<String> filesToDelete,
                              final List<String> filesToEmpty,
                              final int lane,
                              final List<Integer> tiles,
                              final int[] cycles,
                              final String readStructure,
                              final int expectedNumErrors,
                              final List<Integer> filterTiles,
                              final boolean makeFakeFiles) {
        makeFiles(formats, lane, tiles, cycles);
        IlluminaFileUtilTest.deleteRelativeFiles(intensityDir, filesToDelete);
        IlluminaFileUtilTest.emptyRelativeFiles(intensityDir, filesToEmpty);
        writeTileMetricsOutFile(makeMap(makeList(lane - 1, lane + 1, lane), makeList(makeList(1, 2, 3), tiles, tiles)));

        final List<String> args = makeCheckerArgs(basecallDir, lane, readStructure, dataTypes, filterTiles, makeFakeFiles, false);
        Assert.assertEquals(runCommandLine(args), expectedNumErrors);

        // if we previously faked files, then a subsequent run should return no errors!
        if (makeFakeFiles) Assert.assertEquals(runCommandLine(args), 0);
    }

    public void writeFileOfSize(final File file, final int size) {
        try {
            final BufferedWriter writer = new BufferedWriter(new FileWriter(file));
            for (int i = 0; i < size; i++) {
                final int toWrite = Math.min(1000, size);
                final char[] writeBuffer = new char[toWrite];
                for (int j = 0; j < writeBuffer.length; j++) {
                    writeBuffer[j] = (char) (Math.random() * 150);
                }

                writer.write(writeBuffer);
            }
            writer.flush();
            writer.close();
        } catch (final Exception exc) {
            throw new RuntimeException(exc);
        }
    }

    @Test
    public void differentSizedBclTest() {
        final int lane = 5;
        final List<Integer> tiles = makeList(1, 2, 3, 4);
        final int[] cycles = IlluminaFileUtilTest.cycleRange(1, 50);
        final IlluminaDataType[] dataTypes = new IlluminaDataType[]{BaseCalls, IlluminaDataType.QualityScores};

        makeFiles(new SupportedIlluminaFormat[]{Bcl, Filter}, lane, tiles, cycles);
        writeTileMetricsOutFile(makeMap(makeList(lane - 1, lane + 1, lane),
                makeList(makeList(1, 2, 3), tiles, tiles)));

        final File cycleDir = new File(basecallDir, "L005/C9.1");
        writeFileOfSize(new File(cycleDir, "s_5_3.bcl"), 222);

        final List<String> args = makeCheckerArgs(basecallDir, lane, "50T", dataTypes, new ArrayList<Integer>(), false, false);
        Assert.assertEquals(runCommandLine(args), 1);
    }

    @Test(expectedExceptions = SAMException.class)
    public void basedirDoesntExistTest() {
        final List<String> args = makeCheckerArgs(new File("a_made_up_file/in_some_weird_location"), 1, "76T76T",
                new IlluminaDataType[]{IlluminaDataType.Position}, new ArrayList<Integer>(), false, false);
        runCommandLine(args);
    }

    @Test
    public void symlinkLocsTest() {
        final List<Integer> tileList = makeList(1101, 1102, 1103, 2101, 2102, 2103);
        final int lane = 5;
        makeFiles(new SupportedIlluminaFormat[]{Bcl}, lane, tileList, IlluminaFileUtilTest.cycleRange(1, 50));
        List<String> args = makeCheckerArgs(basecallDir, lane, "50T", new IlluminaDataType[]{Position}, new ArrayList<Integer>(), false, true);
        writeTileMetricsOutFile(makeMap(makeList(lane), makeList(tileList)));

        createSingleLocsFile();
        final File intensityLaneDir = new File(intensityDir, IlluminaFileUtil.longLaneStr(lane));
        intensityLaneDir.mkdirs();
        Assert.assertEquals(runCommandLine(args), 0);

        //now that we have created the loc files lets test to make sure they are there
        args = makeCheckerArgs(basecallDir, lane, "50T", new IlluminaDataType[]{IlluminaDataType.Position}, new ArrayList<Integer>(), false, true);
        Assert.assertEquals(runCommandLine(args), 0);
    }

    private void createSingleLocsFile() {
        try {
            final File singleLocsFile = new File(intensityDir, "s.locs");
            final FileWriter writer = new FileWriter(singleLocsFile);
            writer.write("This is a test string.");
            writer.close();
        } catch (final IOException e) {
            e.printStackTrace();
        }
    }
}
