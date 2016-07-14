package org.broadinstitute.hellbender.tools.exome;

import org.apache.commons.io.FileUtils;
import org.apache.commons.io.FilenameUtils;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.ExomeStandardArgumentDefinitions;
import org.broadinstitute.hellbender.tools.exome.acsconversion.ACSModeledSegment;
import org.broadinstitute.hellbender.tools.exome.acsconversion.ACSModeledSegmentUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;


public class CallCNLoHAndSplitsIntegrationTest extends CommandLineProgramTest {
    private static final String ACNV_TEST_SUB_DIR = publicTestDir + "org/broadinstitute/hellbender/tools/exome";
    private static final String TEST_SUB_DIR = publicTestDir + "org/broadinstitute/hellbender/tools/exome/cnlohcaller";

    // Typically, this next file would be the output of a tangent normalization run.
    private static final File COVERAGES_FILE = new File(ACNV_TEST_SUB_DIR, "coverages-for-allelic-integration.tsv");
    private static final File TUMOR_ALLELIC_COUNTS_BASIC_VERBOSITY_FILE = new File(ACNV_TEST_SUB_DIR, "snps-basic.tsv");
    private static final File TUMOR_ALLELIC_COUNTS_FULL_VERBOSITY_FILE = new File(ACNV_TEST_SUB_DIR, "snps-full.tsv");
    private static final File SEGMENT_FILE = new File(TEST_SUB_DIR, "cell_line_small-sim-final.seg");

    @DataProvider(name = "dataCallCNLoHAndSplits")
    public Object[][] dataCallCNLoHAndSplits() {
        return new Object[][]{
                {TUMOR_ALLELIC_COUNTS_BASIC_VERBOSITY_FILE},
                {TUMOR_ALLELIC_COUNTS_FULL_VERBOSITY_FILE}
        };
    }

    @Test(dataProvider = "dataCallCNLoHAndSplits")
    public void testRun(final File tumorAllelicCountsFile) {
        final File outputDir = createTempDir("cnLoH_OutputDir_");
        final List<String> arguments = new ArrayList<>();
        arguments.add("-" + ExomeStandardArgumentDefinitions.TUMOR_ALLELIC_COUNTS_FILE_SHORT_NAME);
        arguments.add(tumorAllelicCountsFile.toString());
        arguments.add("-" + ExomeStandardArgumentDefinitions.TANGENT_NORMALIZED_COUNTS_FILE_SHORT_NAME);
        arguments.add(COVERAGES_FILE.toString());
        arguments.add("-" + ExomeStandardArgumentDefinitions.SEGMENT_FILE_SHORT_NAME);
        arguments.add(SEGMENT_FILE.toString());
        arguments.add("--verbosity");
        arguments.add("INFO");
        arguments.add("-" + CallCNLoHAndSplits.OUTPUT_DIR_SHORT_NAME);
        arguments.add(outputDir.getAbsolutePath());
        arguments.add("-" + CallCNLoHAndSplits.NUM_ITERATIONS_SHORT_NAME);
        arguments.add("3");

        runCommandLine(arguments);

        final String gatkCnvFilename = createFilename(outputDir, CallCNLoHAndSplits.GATK_SEG_FILE_TAG + ".seg", SEGMENT_FILE);
        Assert.assertTrue(new File(gatkCnvFilename).exists());

        final String acsFilename = createFilename(outputDir, CallCNLoHAndSplits.CGA_ACS_SEG_FILE_TAG + ".seg", SEGMENT_FILE);
        Assert.assertTrue(new File(acsFilename).exists());

        final String cnlohFilename = createFilename(outputDir, CallCNLoHAndSplits.CNLOH_BALANCED_SEG_FILE_TAG + ".seg", SEGMENT_FILE);
        Assert.assertTrue(new File(cnlohFilename).exists());

        final String titanTNFilename = createFilename(outputDir, CallCNLoHAndSplits.TITAN_TN_FILE_TAG + ".tsv", SEGMENT_FILE);
        Assert.assertTrue(new File(titanTNFilename).exists());

        final String titanHetFilename = createFilename(outputDir, CallCNLoHAndSplits.TITAN_HET_FILE_TAG + ".tsv", SEGMENT_FILE);
        Assert.assertTrue(new File(titanHetFilename).exists());

        // This is being done to make sure no exception is thrown
        // GATK CNV segs are not log2'd.  ACNV seg means are log2'd
        final List<ModeledSegment> modeledSegments = SegmentUtils.readModeledSegmentsFromSegmentFile(new File(gatkCnvFilename));
        Assert.assertTrue(modeledSegments.size() > 0);

        final int totalTargetCount = modeledSegments.stream().mapToInt(s -> (int)s.getTargetCount()).sum();

        final List<ACSModeledSegment> acsModeledSegments = ACSModeledSegmentUtils.readACSFile(new File(acsFilename));
        Assert.assertTrue(acsModeledSegments.size() > 0);
        final int totalACSTargetCount = acsModeledSegments.stream().mapToInt(s -> (int)s.getTargetCount()).sum();
        Assert.assertEquals(totalACSTargetCount, totalTargetCount);

        final List<ACNVModeledSegment> acnvModeledSegments = SegmentUtils.readACNVModeledSegmentFile(SEGMENT_FILE);
        Assert.assertEquals(acnvModeledSegments.size(), acsModeledSegments.size());

        //TITAN het file cannot be written if original pulldown only has basic verbosity (i.e., is missing ref/alt nucleotide information)
        //in this case, check that a blank file is output
        if (tumorAllelicCountsFile.equals(TUMOR_ALLELIC_COUNTS_BASIC_VERBOSITY_FILE)) {
            try {
                final List<String> titanHetOutputLines = FileUtils.readLines(new File(titanHetFilename));
                Assert.assertEquals(titanHetOutputLines.size(), 0);
            } catch (final IOException e) {
                Assert.fail("Could not read file " + titanHetFilename);
            }
        }
    }

    private String createFilename(final File outputDir, final String extension, final File baseFile) {
        return outputDir.getAbsolutePath() + "/" +
                FilenameUtils.removeExtension(baseFile.getAbsoluteFile().getName()) + "." +
                extension;
    }
}
