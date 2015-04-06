package org.broadinstitute.hellbender.tools.picard.illumina;

import htsjdk.samtools.util.BufferedLineReader;
import htsjdk.samtools.util.LineReader;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileInputStream;
import java.io.PrintWriter;
import java.util.*;

import static htsjdk.samtools.util.StringUtil.join;
import static java.util.Arrays.copyOfRange;
import static org.broadinstitute.hellbender.tools.picard.illumina.IlluminaBasecallsConverter.TILE_NUMBER_COMPARATOR;
import static org.broadinstitute.hellbender.utils.read.SamAssertionUtils.assertSamsEqual;
import static org.testng.Assert.assertTrue;

/**
 * Run IlluminaBasecallsToSam in various barcode & non-barcode modes
 *
 * @author alecw@broadinstitute.org
 */
public class IlluminaBasecallsToSamTest extends CommandLineProgramTest {

    private static final File BASECALLS_DIR = new File(getTestDataDir(), "picard/illumina/25T8B25T/Data/Intensities/BaseCalls");
    private static final File DUAL_BASECALLS_DIR = new File(getTestDataDir(), "picard/illumina/25T8B8B25T/Data/Intensities/BaseCalls");
    private static final File TEST_DATA_DIR = new File(getTestDataDir(), "picard/illumina/25T8B25T/sams");
    private static final File DUAL_TEST_DATA_DIR = new File(getTestDataDir(), "picard/illumina/25T8B8B25T/sams");

    public String getCommandLineProgramName() {
        return IlluminaBasecallsToSam.class.getSimpleName();
    }

    @Test
    public void testTileNumberComparator() {
        assertTrue(TILE_NUMBER_COMPARATOR.compare(100, 10) < 0, "");
        assertTrue(TILE_NUMBER_COMPARATOR.compare(20, 200) > 0, "");
        assertTrue(TILE_NUMBER_COMPARATOR.compare(10, 10) == 0, "");
    }


    @Test
    public void testNonBarcoded() throws Exception {
        runStandardTest(1, "nonBarcode.", "non_barcoded.params", 1, "25S8S25T", BASECALLS_DIR, TEST_DATA_DIR);
    }

    @Test
    public void testMultiplexed() throws Exception {
        runStandardTest(1, "multiplexedBarcode.", "barcode.params", 1, "25T8B25T", BASECALLS_DIR, TEST_DATA_DIR);
    }

    //Same as testMultiplexed except we use BARCODE_1 instead of BARCODE
    @Test
    public void testMultiplexedWithAlternateBarcodeName() throws Exception {
        runStandardTest(1, "singleBarcodeAltName.", "multiplexed_positive_rgtags.params", 1, "25T8B25T", BASECALLS_DIR, TEST_DATA_DIR);
    }

    @Test
    public void testDualBarcodes() throws Exception {
        runStandardTest(1, "dualBarcode.", "barcode_double.params", 2, "25T8B8B25T", DUAL_BASECALLS_DIR, DUAL_TEST_DATA_DIR);
    }

    /**
     * Ensures that a run missing a barcode from the parameters file throws an error.
     * <p>
     * TODO: This testcase isn't broken, but can spawn an issue with FileChannelJDKBugWorkAround since it expects
     * an exception to be thrown.
     */
    @Test(groups = {"broken"})
    public void testCorruptDataReturnCode() throws Exception {
        boolean exceptionThrown = false;
        try {
            runStandardTest(9, "dualBarcode.", "negative_test.params", 2, "30T8B8B", BASECALLS_DIR, TEST_DATA_DIR);
        } catch (Throwable e) {
            exceptionThrown = true;
        } finally {
            assertTrue(exceptionThrown);
        }
    }

    /**
     * This test utility takes a libraryParamsFile and generates output sam files through IlluminaBasecallsToSam to compare against
     * preloaded test data
     *
     * @param jobName
     * @param libraryParamsFile
     * @param concatNColumnFields
     * @param readStructure
     * @throws Exception
     */
    private void runStandardTest(final int lane, final String jobName, final String libraryParamsFile,
                                 final int concatNColumnFields, final String readStructure,
                                 final File baseCallsDir, final File testDataDir) throws Exception {
        final File outputDir = File.createTempFile(jobName, ".dir");
        outputDir.delete();
        outputDir.mkdir();
        outputDir.deleteOnExit();
        // Create barcode.params with output files in the temp directory
        final File libraryParams = new File(outputDir, libraryParamsFile);
        libraryParams.deleteOnExit();
        final List<File> samFiles = new ArrayList<File>();
        final LineReader reader = new BufferedLineReader(new FileInputStream(new File(testDataDir, libraryParamsFile)));
        final PrintWriter writer = new PrintWriter(libraryParams);
        final String header = reader.readLine();
        writer.println(header + "\tOUTPUT");
        String line = reader.readLine();
        while (line != null) {
            final String[] fields = line.split("\t");
            final File outputSam = new File(outputDir, join("", copyOfRange(fields, 0, concatNColumnFields)) + ".sam");
            outputSam.deleteOnExit();
            samFiles.add(outputSam);
            writer.println(line + "\t" + outputSam);
            line = reader.readLine();
        }
        writer.close();
        reader.close();

        runCommandLine(new String[]{
                "--BASECALLS_DIR", baseCallsDir.getAbsolutePath(),
                "--LANE", Integer.toString(lane),
                "--RUN_BARCODE", "HiMom",
                "--READ_STRUCTURE", readStructure,
                "--LIBRARY_PARAMS", libraryParams.getAbsolutePath()
        });

        for (final File outputSam : samFiles) {
            assertSamsEqual(outputSam, new File(testDataDir, outputSam.getName()));
        }
    }
}
