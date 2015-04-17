package org.broadinstitute.hellbender.tools.picard.illumina;

import htsjdk.samtools.util.BufferedLineReader;
import htsjdk.samtools.util.LineReader;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.tools.picard.illumina.parser.ReadStructure;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileInputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;

import static htsjdk.samtools.util.IOUtil.assertFilesEqual;
import static htsjdk.samtools.util.StringUtil.join;
import static htsjdk.samtools.util.TestUtil.recursiveDelete;
import static java.util.Arrays.copyOfRange;
import static org.broadinstitute.hellbender.tools.picard.illumina.IlluminaBasecallsToFastq.ReadNameFormat.ILLUMINA;

public class IlluminaBasecallsToFastqTest extends CommandLineProgramTest {

    private static final File BASECALLS_DIR = new File(getTestDataDir(), "picard/illumina/25T8B25T/Data/Intensities/BaseCalls");
    private static final File DUAL_BASECALLS_DIR = new File(getTestDataDir(), "/picard/illumina/25T8B8B25T/Data/Intensities/BaseCalls");
    private static final File TEST_DATA_DIR = new File(getTestDataDir(), "picard/illumina/25T8B25T/fastq");
    private static final File DUAL_TEST_DATA_DIR = new File(getTestDataDir(), "picard/illumina/25T8B8B25T/fastq");

    public String getCommandLineProgramName() {
        return IlluminaBasecallsToFastq.class.getSimpleName();
    }

    @Test(enabled = false, description = "bug https://github.com/broadinstitute/hellbender/issues/364")
    public void testNonBarcoded() throws Exception {
        final String suffix = ".1.fastq";
        final File outputFastq1 = File.createTempFile("nonBarcoded.", suffix);
        outputFastq1.deleteOnExit();
        final String outputPrefix = outputFastq1.getAbsolutePath().substring(0, outputFastq1.getAbsolutePath().length() - suffix.length());
        final File outputFastq2 = new File(outputPrefix + ".2.fastq");
        outputFastq2.deleteOnExit();
        final int lane = 1;
        runCommandLine(new String[]{
                "--BASECALLS_DIR", BASECALLS_DIR.getAbsolutePath(),
                "--LANE", Integer.toString(lane),
                "--READ_STRUCTURE", "25T8B25T",
                "--OUTPUT_PREFIX", outputPrefix,
                "--RUN_BARCODE", "HiMom",
                "--MACHINE_NAME", "machine1",
                "--FLOWCELL_BARCODE", "abcdeACXX"
        });
        assertFilesEqual(outputFastq1, new File(TEST_DATA_DIR, "nonBarcoded.1.fastq"));
        assertFilesEqual(outputFastq2, new File(TEST_DATA_DIR, "nonBarcoded.2.fastq"));
    }

    @Test(enabled = false, description = "bug https://github.com/broadinstitute/hellbender/issues/364")
    public void testMultiplexWithIlluminaReadNameHeaders() throws Exception {
        final File outputDir = File.createTempFile("testMultiplexRH.", ".dir");
        try {
            outputDir.delete();
            outputDir.mkdir();
            outputDir.deleteOnExit();

            final String filePrefix = "testMultiplexRH";
            final File outputPrefix = new File(outputDir, filePrefix);

            runCommandLine(new String[]{
                    "--BASECALLS_DIR", BASECALLS_DIR.getAbsolutePath(),
                    "--LANE", "1",
                    "--RUN_BARCODE", "HiMom",
                    "--READ_STRUCTURE", "25T8B25T",
                    "--OUTPUT_PREFIX", outputPrefix.getAbsolutePath(),
                    "--MACHINE_NAME", "machine1",
                    "--FLOWCELL_BARCODE", "abcdeACXX",
                    "--READ_NAME_FORMAT", ILLUMINA.toString()
            });

            final String[] filenames = new String[]{
                    filePrefix + ".1.fastq",
                    filePrefix + ".barcode_1.fastq"
            };
            for (final String filename : filenames) {
                assertFilesEqual(new File(outputDir, filename), new File(TEST_DATA_DIR, filename));
            }

        } finally {
            recursiveDelete(outputDir);
        }
    }

    @Test(enabled = false, description = "bug https://github.com/broadinstitute/hellbender/issues/364")
    public void testDeMultiplexed() throws Exception {
        runStandardTest(1, "multiplexedBarcode.", "mp_barcode.params", 1, "25T8B25T", BASECALLS_DIR, TEST_DATA_DIR);
    }

    @Test(enabled = false, description = "bug https://github.com/broadinstitute/hellbender/issues/364")
    public void testDualBarcodes() throws Exception {
        runStandardTest(1, "dualBarcode.", "barcode_double.params", 2, "25T8B8B25T", DUAL_BASECALLS_DIR, DUAL_TEST_DATA_DIR);
    }

    /**
     * This test utility takes a libraryParamsFile and generates output sam files through IlluminaBasecallsToFastq to compare against
     * preloaded test data
     *
     * @param jobName
     * @param libraryParamsFile
     * @param concatNColumnFields
     * @param readStructureString
     * @throws Exception
     */
    private void runStandardTest(final int lane, final String jobName, final String libraryParamsFile,
                                 final int concatNColumnFields, final String readStructureString, final File baseCallsDir,
                                 final File testDataDir) throws Exception {
        final File outputDir = File.createTempFile(jobName, ".dir");
        try {
            outputDir.delete();
            outputDir.mkdir();
            outputDir.deleteOnExit();
            // Create barcode.params with output files in the temp directory
            final File libraryParams = new File(outputDir, libraryParamsFile);
            libraryParams.deleteOnExit();
            final List<File> outputPrefixes = new ArrayList<File>();
            final LineReader reader = new BufferedLineReader(new FileInputStream(new File(testDataDir, libraryParamsFile)));
            final PrintWriter writer = new PrintWriter(libraryParams);
            final String header = reader.readLine();
            writer.println(header + "\tOUTPUT_PREFIX");
            String line = reader.readLine();
            while (line != null) {
                final String[] fields = line.split("\t");
                final File outputPrefix = new File(outputDir, join("", copyOfRange(fields, 0, concatNColumnFields)));
                outputPrefixes.add(outputPrefix);
                writer.println(line + "\t" + outputPrefix);
                line = reader.readLine();
            }
            writer.close();
            reader.close();

            runCommandLine(new String[]{
                    "--BASECALLS_DIR", baseCallsDir.getAbsolutePath(),
                    "--LANE", Integer.toString(lane),
                    "--RUN_BARCODE", "HiMom",
                    "--READ_STRUCTURE", readStructureString,
                    "--MULTIPLEX_PARAMS", libraryParams.getAbsolutePath(),
                    "--MACHINE_NAME", "machine1",
                    "--FLOWCELL_BARCODE", "abcdeACXX"
            });

            final ReadStructure readStructure = new ReadStructure(readStructureString);
            for (final File outputSam : outputPrefixes) {
                for (int i = 1; i <= readStructure.templates.length(); ++i) {
                    final String filename = outputSam.getName() + "." + i + ".fastq";
                    assertFilesEqual(new File(outputSam.getParentFile(), filename), new File(testDataDir, filename));
                }
                for (int i = 1; i <= readStructure.barcodes.length(); ++i) {
                    final String filename = outputSam.getName() + ".barcode_" + i + ".fastq";
                    assertFilesEqual(new File(outputSam.getParentFile(), filename), new File(testDataDir, filename));
                }
            }
        } finally {
            recursiveDelete(outputDir);
        }
    }
}
