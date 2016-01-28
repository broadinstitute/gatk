package org.broadinstitute.hellbender.tools.exome;

import org.apache.commons.io.FileUtils;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.ExomeStandardArgumentDefinitions;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class AllelicCNVIntegrationTest extends CommandLineProgramTest {
    private static final String TEST_SUB_DIR = publicTestDir + "org/broadinstitute/hellbender/tools/exome/";

    private static final File COVERAGES_FILE = new File(TEST_SUB_DIR
            + "coverages-for-allelic-integration.tsv");
    private static final File SNP_COUNTS_FILE = new File(TEST_SUB_DIR
            + "snps-for-allelic-integration.tsv");
    private static final File SEGMENT_FILE =
            new File(TEST_SUB_DIR + "segments-for-allelic-integration.seg");
    private static final String SAMPLE_NAME = "test";

    @Test
    public void testACNV() {
        final File tempDir = createTempDir("allelic-integration-" + SAMPLE_NAME);
        final String tempDirPath = tempDir.getAbsolutePath();
        final String outputPrefix = tempDirPath + "/" + SAMPLE_NAME;

        final String[] arguments = {
                "--" + ExomeStandardArgumentDefinitions.TUMOR_ALLELIC_COUNTS_FILE_LONG_NAME, SNP_COUNTS_FILE.getAbsolutePath(),
                "--" + ExomeStandardArgumentDefinitions.TANGENT_NORMALIZED_COUNTS_LONG_NAME, COVERAGES_FILE.getAbsolutePath(),
                "--" + ExomeStandardArgumentDefinitions.SEGMENT_FILE_LONG_NAME, SEGMENT_FILE.getAbsolutePath(),
                "--" + AllelicCNV.OUTPUT_PREFIX_LONG_NAME, outputPrefix,
                "--" + ExomeStandardArgumentDefinitions.SAMPLE_LONG_NAME, SAMPLE_NAME,
                "--" + AllelicCNV.NUM_SAMPLES_COPY_RATIO_LONG_NAME, "25",
                "--" + AllelicCNV.NUM_BURN_IN_COPY_RATIO_LONG_NAME, "10",
                "--" + AllelicCNV.NUM_SAMPLES_ALLELE_FRACTION_LONG_NAME, "25",
                "--" + AllelicCNV.NUM_BURN_IN_ALLELE_FRACTION_LONG_NAME, "10",
                "--verbosity", "INFO",
        };
        runCommandLine(arguments);

        //only check that files are created, do not check for correctness of results
        final File snpSegmentsFile = new File(outputPrefix + "-" + AllelicCNV.SNP_MAF_SEG_FILE_TAG + ".seg");
        final File unionedSegmentsFile = new File(outputPrefix + "-" + AllelicCNV.UNION_SEG_FILE_TAG + ".seg");
        final File noSmallSegmentsFile = new File(outputPrefix + "-" + AllelicCNV.SMALL_MERGED_SEG_FILE_TAG + ".seg");
        final File initialSimilarSegmentsFile = new File(outputPrefix + "-" + ACNVModeller.INITIAL_SEG_FILE_TAG + ".seg");
        final File finalSimilarSegmentsFile = new File(outputPrefix + "-" + ACNVModeller.FINAL_SEG_FILE_TAG + ".seg");

        final List<File> outputFileList = new ArrayList<>();
        outputFileList.add(snpSegmentsFile);
        outputFileList.add(unionedSegmentsFile);
        outputFileList.add(noSmallSegmentsFile);
        outputFileList.add(initialSimilarSegmentsFile);
        outputFileList.add(finalSimilarSegmentsFile);

        for (final File outputFile: outputFileList) {

            // Check that all files are files with a size greater than 0.
            Assert.assertTrue(outputFile.isFile());
            Assert.assertTrue(outputFile.length() > 0);
            try {

                //Check that all files have:
                //  - at least two lines and all either start with "#" or contain at least one "\t"
                //  - at least two lines with tab (header + 1 segment)
                final List<String> outputLines = FileUtils.readLines(outputFile);
                Assert.assertTrue(outputLines.size() > 2);
                Assert.assertEquals(outputLines.stream().filter(l -> l.contains("\t") || l.startsWith("#")).count(), outputLines.size());
                Assert.assertTrue(outputLines.stream().filter(l -> l.split("\t").length > 2).count() > 2, "File: " + outputFile + " does not seem to have at least one segment and a header.");

            } catch (final IOException ioe) {
                Assert.fail("Could not read file: " + outputFile, ioe);
            }
        }
    }
}