package org.broadinstitute.hellbender.tools.exome;

import org.apache.commons.io.filefilter.WildcardFileFilter;
import org.apache.commons.lang3.ArrayUtils;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.tools.pon.allelic.AllelicPanelOfNormals;
import org.broadinstitute.hellbender.tools.pon.allelic.AllelicPoNTestUtils;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileFilter;
import java.util.Arrays;
import java.util.Collection;
import java.util.stream.Stream;

import static org.broadinstitute.hellbender.tools.exome.CreateAllelicPanelOfNormals.SITE_FREQUENCY_THRESHOLD_LONG_NAME;
import static org.broadinstitute.hellbender.tools.exome.CreateAllelicPanelOfNormals.TSV_OUTPUT_FILE_LONG_NAME;

/**
 * Integration tests for {@link CreateAllelicPanelOfNormals}.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public class CreateAllelicPanelOfNormalsIntegrationTest extends CommandLineProgramTest {
    private static final String TEST_SUB_DIR = publicTestDir + "org/broadinstitute/hellbender/tools/exome/";

    private static final FileFilter PULLDOWN_FILTER = new WildcardFileFilter("allelic-pon-test-pulldown-*tsv");
    private static final File[] PULLDOWN_FILES = new File(TEST_SUB_DIR).listFiles(PULLDOWN_FILTER);
    private static final String[] PULLDOWN_FILE_ARGUMENTS = Stream.of(PULLDOWN_FILES)
            .map(f -> Arrays.asList("--" + StandardArgumentDefinitions.INPUT_LONG_NAME, f.getAbsolutePath()))
            .flatMap(Collection::stream)
            .toArray(size -> new String[2 * PULLDOWN_FILES.length]);

    private static final File PON_EXPECTED_SITE_FREQUENCY_50_TSV = new File(TEST_SUB_DIR, "allelic-pon-test-pon-freq-50.tsv");
    private static final File PON_EXPECTED_SITE_FREQUENCY_75_TSV = new File(TEST_SUB_DIR, "allelic-pon-test-pon-freq-75.tsv");

    private static final File PON_EXPECTED_SITE_FREQUENCY_50_HDF5 = new File(TEST_SUB_DIR, "allelic-pon-test-pon-freq-50.pon");
    private static final File PON_EXPECTED_SITE_FREQUENCY_75_HDF5 = new File(TEST_SUB_DIR, "allelic-pon-test-pon-freq-75.pon");

    private static final String TEST_HDF5_OUTPUT_FILE_NAME = "create-allelic-pon-test.pon";

    @DataProvider(name = "dataArgumentValidation")
    public Object[][] dataArgumentValidation() {
        return new Object[][]{
                {new String[]{"--" + SITE_FREQUENCY_THRESHOLD_LONG_NAME, "0"}},
                {new String[]{"--" + SITE_FREQUENCY_THRESHOLD_LONG_NAME, "10"}},
                {new String[]{"--" + SITE_FREQUENCY_THRESHOLD_LONG_NAME, "-1"}}
        };
    }

    @Test(dataProvider = "dataArgumentValidation", expectedExceptions = IllegalArgumentException.class)
    public void testArgumentValidation(final String[] testArguments) {
        final String[] commonArguments = ArrayUtils.addAll(PULLDOWN_FILE_ARGUMENTS,
                "--" + StandardArgumentDefinitions.OUTPUT_LONG_NAME, TEST_HDF5_OUTPUT_FILE_NAME);
        final String[] arguments = ArrayUtils.addAll(commonArguments, testArguments);
        runCommandLine(arguments);
    }

    @DataProvider(name = "dataSiteFrequency")
    public Object[][] dataSiteFrequency() {
        return new Object[][]{
                {new String[]{"--" + SITE_FREQUENCY_THRESHOLD_LONG_NAME, "0.5"}, PON_EXPECTED_SITE_FREQUENCY_50_HDF5, PON_EXPECTED_SITE_FREQUENCY_50_TSV},
                {new String[]{"--" + SITE_FREQUENCY_THRESHOLD_LONG_NAME, "0.75"}, PON_EXPECTED_SITE_FREQUENCY_75_HDF5, PON_EXPECTED_SITE_FREQUENCY_75_TSV}
        };
    }

    @Test(dataProvider = "dataSiteFrequency")
    public void testHDF5OutputOnly(final String[] testArguments, final File expectedHDF5File, final File expectedTSVFile) {
        final File allelicPoNHDF5File = createTempFile("create-allelic-pon-test", ".pon");
        allelicPoNHDF5File.delete();
        final String[] commonArguments = ArrayUtils.addAll(PULLDOWN_FILE_ARGUMENTS,
                "--" + StandardArgumentDefinitions.OUTPUT_LONG_NAME, allelicPoNHDF5File.getAbsolutePath());
        final String[] arguments = ArrayUtils.addAll(commonArguments, testArguments);

        runCommandLine(arguments);
        final AllelicPanelOfNormals resultHDF5 = AllelicPanelOfNormals.read(allelicPoNHDF5File);
        final AllelicPanelOfNormals expectedHDF5 = AllelicPanelOfNormals.read(expectedHDF5File);
        AllelicPoNTestUtils.assertAllelicPoNsEqual(resultHDF5, expectedHDF5);

        //check overwrite
        runCommandLine(arguments);
        final AllelicPanelOfNormals resultHDF5Overwrite = AllelicPanelOfNormals.read(allelicPoNHDF5File);
        AllelicPoNTestUtils.assertAllelicPoNsEqual(resultHDF5Overwrite, expectedHDF5);
    }

    @Test(dataProvider = "dataSiteFrequency", groups = "createsTempFiles")
    public void testHDF5AndTSVOutput(final String[] testArguments, final File expectedHDF5File, final File expectedTSVFile) {
        final File allelicPoNHDF5File = createTempFile("create-allelic-pon-test", ".pon");
        allelicPoNHDF5File.delete();
        final File allelicPoNTSVFile = createTempFile("create-allelic-pon-test", ".tsv");
        allelicPoNTSVFile.delete();

        final String[] commonArguments = ArrayUtils.addAll(PULLDOWN_FILE_ARGUMENTS,
                "--" + StandardArgumentDefinitions.OUTPUT_LONG_NAME, allelicPoNHDF5File.getAbsolutePath(),
                "--" + TSV_OUTPUT_FILE_LONG_NAME, allelicPoNTSVFile.getAbsolutePath());
        final String[] arguments = ArrayUtils.addAll(commonArguments, testArguments);

        runCommandLine(arguments);
        final AllelicPanelOfNormals resultHDF5 = AllelicPanelOfNormals.read(allelicPoNHDF5File);
        final AllelicPanelOfNormals expectedHDF5 = AllelicPanelOfNormals.read(expectedHDF5File);
        AllelicPoNTestUtils.assertAllelicPoNsEqual(resultHDF5, expectedHDF5);
        final AllelicPanelOfNormals resultTSV = AllelicPanelOfNormals.read(allelicPoNTSVFile);
        final AllelicPanelOfNormals expectedTSV = AllelicPanelOfNormals.read(expectedTSVFile);
        AllelicPoNTestUtils.assertAllelicPoNsEqual(resultTSV, expectedTSV);

        //check overwrite
        runCommandLine(arguments);
        final AllelicPanelOfNormals resultHDF5Overwrite = AllelicPanelOfNormals.read(allelicPoNHDF5File);
        AllelicPoNTestUtils.assertAllelicPoNsEqual(resultHDF5Overwrite, expectedHDF5);
        final AllelicPanelOfNormals resultTSVOverwrite = AllelicPanelOfNormals.read(allelicPoNTSVFile);
        AllelicPoNTestUtils.assertAllelicPoNsEqual(resultTSVOverwrite, expectedTSV);
    }
}