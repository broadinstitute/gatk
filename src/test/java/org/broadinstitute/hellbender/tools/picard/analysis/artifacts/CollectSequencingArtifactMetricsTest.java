package org.broadinstitute.hellbender.tools.picard.analysis.artifacts;

import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.IOUtil;
import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.testng.Assert;
import org.testng.annotations.AfterTest;
import org.testng.annotations.BeforeTest;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

public final class CollectSequencingArtifactMetricsTest extends CommandLineProgramTest {

    private static final File TEST_DIR = new File(getTestDataDir(), "picard/analysis/artifacts/CollectSequencingArtifactMetrics");
    private static final File REFERENCE = new File(TEST_DIR, "test.fasta");
    private static final File TEST_SAM = new File(TEST_DIR, "test.sam");
    private static final File DB_SNP = new File(TEST_DIR, "test.dbsnp.vcf");
    private static final File INTERVALS = new File(TEST_DIR, "test.interval_list");
    private static final File TEST_CASES = new File(TEST_DIR, "ExpectedMetricsOutput");

    private File globalTempOutputDir;

    @Override
    public String getTestedClassName() {
        return CollectSequencingArtifactMetrics.class.getSimpleName();
    }

    @BeforeTest
    public void setUp() throws IOException {
        globalTempOutputDir = IOUtil.createTempDir("artifactMetrics.", ".tmp");
    }

    @AfterTest
    public void tearDown() throws IOException {
        IOUtil.deleteDirectoryTree(globalTempOutputDir);
    }

    /**
     * Run the CLP using standard arguments and maybe some additional ones, then compare to the expected results.
     *
     * @param testCase name of test case (should match one of the file sets in {@code TEST_CASES}
     * @param extraArgs optional arguments, spaced-delimited, as they would be passed to the command line
     */
    @Test(dataProvider = "data")
    public void runAnalysis(final String testCase, final String extraArgs) throws IOException {
        final File actual = new File(globalTempOutputDir, testCase);
        final File expected = new File(TEST_CASES, testCase);

        final StringBuilder args = new StringBuilder();
        args.append("--input " + TEST_SAM.getAbsolutePath());
        args.append(" --output " + actual.getAbsolutePath());
        args.append(" --reference " + REFERENCE.getAbsolutePath());

        if (!extraArgs.isEmpty()) args.append(" " + extraArgs);

        runCommandLine(StringUtils.split(args.toString()));
        assertAllFilesEqual(expected, actual);
    }

    private void assertAllFilesEqual(final File expectedBase, final File actualBase) {
        boolean equal = areMetricsEqual(expectedBase, actualBase, SequencingArtifactMetrics.PRE_ADAPTER_SUMMARY_EXT);
        equal = equal && areMetricsEqual(expectedBase, actualBase, SequencingArtifactMetrics.PRE_ADAPTER_DETAILS_EXT);
        equal = equal && areMetricsEqual(expectedBase, actualBase, SequencingArtifactMetrics.BAIT_BIAS_SUMMARY_EXT);
        equal = equal && areMetricsEqual(expectedBase, actualBase, SequencingArtifactMetrics.BAIT_BIAS_DETAILS_EXT);
        Assert.assertTrue(equal);
    }

    private boolean areMetricsEqual(final File expectedBase, final File actualBase, final String extension) {
        return MetricsFile.areMetricsEqual(new File(expectedBase + extension), new File(actualBase + extension));
    }

    @DataProvider(name = "data")
    public Object[][] testData() {
        // default CONTEXT_SIZE is 0 to cut down on file sizes
        return new Object[][]{
                {"with_context",   "--MIN_INS 30 --MAX_INS 30 --CONTEXT_SIZE 1"},
                {"with_dbsnp",     "--MIN_INS 30 --MAX_INS 30 --CONTEXT_SIZE 0 --DB_SNP " + DB_SNP},
                {"with_intervals", "--MIN_INS 30 --MAX_INS 30 --CONTEXT_SIZE 0 --INTERVALS " + INTERVALS},
                {"no_bq_cutoff",   "--MIN_INS 30 --MAX_INS 30 --CONTEXT_SIZE 0 --MINIMUM_QUALITY_SCORE 0"},
                {"no_mq_cutoff",   "--MIN_INS 30 --MAX_INS 30 --CONTEXT_SIZE 0 --MINIMUM_MAPPING_QUALITY 0"},
                {"unmapped_mate",  "--MIN_INS 0  --MAX_INS 0  --CONTEXT_SIZE 0 --MINIMUM_MAPPING_QUALITY 0"},
        };
    }
}
