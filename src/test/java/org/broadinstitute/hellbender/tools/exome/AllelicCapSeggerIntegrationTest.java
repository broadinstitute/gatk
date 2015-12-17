package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.ExomeStandardArgumentDefinitions;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;

public class AllelicCapSeggerIntegrationTest extends CommandLineProgramTest {
    private static final String TEST_SUB_DIR = publicTestDir + "org/broadinstitute/hellbender/tools/exome/";

    private static final File COVERAGES_FILE = new File(TEST_SUB_DIR
            + "coverages-for-allelic-integration.tsv");
    private static final File SNP_COUNTS_FILE = new File(TEST_SUB_DIR
            + "snps-for-allelic-integration.tsv");
    private static final File SEGMENT_FILE =
            new File(TEST_SUB_DIR + "segments-for-allelic-integration.seg");
    private static final String SAMPLE_NAME = "test";

    @Test
    public void testAllelicCapSegger() {
        final File tempDir = createTempDir("allelic-integration-" + SAMPLE_NAME);
        final String tempDirPath = tempDir.getAbsolutePath();
        final String outputPrefix = tempDirPath + "/" + SAMPLE_NAME;

        final String[] arguments = {
                "--" + ExomeStandardArgumentDefinitions.TUMOR_ALLELIC_COUNTS_FILE_LONG_NAME, SNP_COUNTS_FILE.getAbsolutePath(),
                "--" + ExomeStandardArgumentDefinitions.TANGENT_NORMALIZED_COUNTS_LONG_NAME, COVERAGES_FILE.getAbsolutePath(),
                "--" + ExomeStandardArgumentDefinitions.SEGMENT_FILE_LONG_NAME, SEGMENT_FILE.getAbsolutePath(),
                "--" + AllelicCapSegger.OUTPUT_PREFIX_LONG_NAME, outputPrefix,
                "--" + ExomeStandardArgumentDefinitions.SAMPLE_LONG_NAME, SAMPLE_NAME,
                "--" + AllelicCapSegger.NUM_SAMPLES_COPY_RATIO_LONG_NAME, "25",
                "--" + AllelicCapSegger.NUM_BURN_IN_COPY_RATIO_LONG_NAME, "10",
                "--" + AllelicCapSegger.NUM_SAMPLES_ALLELE_FRACTION_LONG_NAME, "25",
                "--" + AllelicCapSegger.NUM_BURN_IN_ALLELE_FRACTION_LONG_NAME, "10",
                "--VERBOSITY", "INFO",
        };
        runCommandLine(arguments);

        //only check that files are created, do not check for correctness of results
        final File snpSegmentsFile = new File(outputPrefix + "-" + AllelicCapSegger.SNP_MAF_SEG_FILE_TAG + ".seg");
        final File unionedSegmentsFile = new File(outputPrefix + "-" + AllelicCapSegger.UNION_SEG_FILE_TAG + ".seg");
        final File noSmallSegmentsFile = new File(outputPrefix + "-" + AllelicCapSegger.SMALL_MERGED_SEG_FILE_TAG + ".seg");
        final File initialSimilarSegmentsFile = new File(outputPrefix + "-" + ACSModeller.INITIAL_SEG_FILE_TAG + ".seg");
        final File finalSimilarSegmentsFile = new File(outputPrefix + "-" + ACSModeller.FINAL_SEG_FILE_TAG + ".seg");

        Assert.assertTrue(snpSegmentsFile.isFile());
        Assert.assertTrue(unionedSegmentsFile.isFile());
        Assert.assertTrue(noSmallSegmentsFile.isFile());
        Assert.assertTrue(initialSimilarSegmentsFile.isFile());
        Assert.assertTrue(finalSimilarSegmentsFile.isFile());
    }
}