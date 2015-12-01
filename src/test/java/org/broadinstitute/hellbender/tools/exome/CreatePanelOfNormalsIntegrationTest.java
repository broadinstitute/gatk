package org.broadinstitute.hellbender.tools.exome;

import htsjdk.samtools.util.Log;
import org.apache.commons.math3.linear.RealMatrix;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.LoggingUtils;
import org.broadinstitute.hellbender.utils.hdf5.HDF5File;
import org.broadinstitute.hellbender.utils.hdf5.HDF5PoN;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.LinkedHashSet;
import java.util.List;

/**
 * Integration test for {@link CreatePanelOfNormals}.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public class CreatePanelOfNormalsIntegrationTest extends CommandLineProgramTest {

    private static final File TEST_FILE_DIR = new File("src/test/resources/org/broadinstitute/hellbender/tools/exome");

    private static final File CONTROL_PCOV_FULL_FILE = new File(TEST_FILE_DIR, "create-pon-control-full.pcov");

    private static final File CONTROL_PCOV_SOME_TARGETS_FULL_FILE = new File(TEST_FILE_DIR, "create-pon-control-some-targets-full.pcov");

    private static final File CONTROL_PCOV_TARGET_NAME_ONLY_FILE = new File(TEST_FILE_DIR, "create-pon-control-target-name-only.pcov");

    private static final File CONTROL_PCOV_TARGET_COORDINATE_ONLY_FILE = new File(TEST_FILE_DIR, "create-pon-control-target-coord-only.pcov");

    private static final File ALL_TARGETS_FILE = new File(TEST_FILE_DIR, "create-pon-all-targets.tab");

    private static final File SOME_TARGETS_FILE = new File(TEST_FILE_DIR, "create-pon-some-targets.tab");

    private static final File EXPECTED_ALL_TARGETS_PON = new File(TEST_FILE_DIR, "create-pon-all-targets.pon");

    private static final File EXPECTED_SOME_TARGETS_PON = new File(TEST_FILE_DIR, "create-pon-some-targets.pon");

    @Test(dataProvider="allTargetsHDF5PoNCreationData")
    public void testAllTargetsHDF5PoNCreationNoSpark(final File targetsFile, final File inputFile) {
        final List<String> arguments = new ArrayList<>();
        arguments.add("-" + StandardArgumentDefinitions.INPUT_SHORT_NAME);
        arguments.add(inputFile.toString());
        if (targetsFile != null) {
            arguments.add("-" + TargetArgumentCollection.TARGET_FILE_SHORT_NAME);
            arguments.add(targetsFile.toString());
        }
        final File outputFile = createTempFile("pon-",".hd5");
        arguments.add("-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME);
        arguments.add(outputFile.toString());
        arguments.add("--VERBOSITY");
        arguments.add("INFO");
        arguments.add("-ds");

        runCommandLine(arguments);
        assertEquivalentPoN(outputFile, EXPECTED_ALL_TARGETS_PON);
    }

    @Test(dataProvider="allTargetsHDF5PoNCreationData")
    public void testAllTargetsHDF5PoNCreationSpark(final File targetsFile, final File inputFile) {
        final List<String> arguments = new ArrayList<>();
        arguments.add("-" + StandardArgumentDefinitions.INPUT_SHORT_NAME);
        arguments.add(inputFile.toString());
        if (targetsFile != null) {
            arguments.add("-" + TargetArgumentCollection.TARGET_FILE_SHORT_NAME);
            arguments.add(targetsFile.toString());
        }
        final File outputFile = createTempFile("pon-",".hd5");
        arguments.add("-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME);
        arguments.add(outputFile.toString());
        arguments.add("--VERBOSITY");
        arguments.add("INFO");

        runCommandLine(arguments);
        assertEquivalentPoN(outputFile, EXPECTED_ALL_TARGETS_PON);
    }


    @Test()
    public void testDryRun() {
        final List<String> arguments = new ArrayList<>();
        arguments.add("-" + StandardArgumentDefinitions.INPUT_SHORT_NAME);
        arguments.add(CONTROL_PCOV_FULL_FILE.toString());
        final File outputFile = createTempFile("pon-",".hd5");
        outputFile.delete();
        arguments.add("-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME);
        arguments.add(outputFile.toString());
        arguments.add("-" + CreatePanelOfNormals.DRY_RUN_SHORT_NAME);
        runCommandLine(arguments);
        Assert.assertFalse(outputFile.exists());
    }

    @Test(expectedExceptions = UserException.class)
    public void testCoordinatesOnly() {
        final List<String> arguments = new ArrayList<>();
        arguments.add("-" + StandardArgumentDefinitions.INPUT_SHORT_NAME);
        arguments.add(CONTROL_PCOV_TARGET_COORDINATE_ONLY_FILE.toString());
        final File outputFile = createTempFile("pon-",".hd5");
        arguments.add("-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME);
        arguments.add(outputFile.toString());
        runCommandLine(arguments);
    }

    @Test(dataProvider="badPercentileMax50Data", expectedExceptions = UserException.BadArgumentValue.class)
    public void testBadCountTruncatePercentile(final double percentile) {
        final List<String> arguments = new ArrayList<>();
        arguments.add("-" + StandardArgumentDefinitions.INPUT_SHORT_NAME);
        arguments.add(CONTROL_PCOV_FULL_FILE.toString());
        final File outputFile = createTempFile("pon-",".hd5");
        arguments.add("-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME);
        arguments.add(outputFile.toString());
        arguments.add("-" + CreatePanelOfNormals.COUNT_TRUNCATE_PERCENTILE_SHORT_NAME);
        arguments.add("" + percentile);
        runCommandLine(arguments);
    }

    @Test(dataProvider="badPercentileMax100Data", expectedExceptions = UserException.BadArgumentValue.class)
    public void testBadMaximumAllowedZerosInColumn(final double percentile) {
        final List<String> arguments = new ArrayList<>();
        arguments.add("-" + StandardArgumentDefinitions.INPUT_SHORT_NAME);
        arguments.add(CONTROL_PCOV_FULL_FILE.toString());
        final File outputFile = createTempFile("pon-",".hd5");
        arguments.add("-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME);
        arguments.add(outputFile.toString());
        arguments.add("-" + CreatePanelOfNormals.MAXIMUM_PERCENT_ZEROS_IN_COLUMN_SHORT_NAME);
        arguments.add("" + percentile);
        runCommandLine(arguments);
    }

    @Test(dataProvider="badPercentileMax100Data", expectedExceptions = UserException.BadArgumentValue.class)
    public void testBadMaximumAllowedZerosInTarget(final double percentile) {
        final List<String> arguments = new ArrayList<>();
        arguments.add("-" + StandardArgumentDefinitions.INPUT_SHORT_NAME);
        arguments.add(CONTROL_PCOV_FULL_FILE.toString());
        final File outputFile = createTempFile("pon-",".hd5");
        arguments.add("-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME);
        arguments.add(outputFile.toString());
        arguments.add("-" + CreatePanelOfNormals.MAXIMUM_PERCENT_ZEROS_IN_TARGET_SHORT_NAME);
        arguments.add("" + percentile);
        runCommandLine(arguments);
    }

    @Test(dataProvider="badPercentileMax100Data", expectedExceptions = UserException.BadArgumentValue.class)
    public void testBadPercentileUsableTargets(final double percentile) {
        final List<String> arguments = new ArrayList<>();
        arguments.add("-" + StandardArgumentDefinitions.INPUT_SHORT_NAME);
        arguments.add(CONTROL_PCOV_FULL_FILE.toString());
        final File outputFile = createTempFile("pon-",".hd5");
        arguments.add("-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME);
        arguments.add(outputFile.toString());
        arguments.add("-" + CreatePanelOfNormals.TARGET_FACTOR_THRESHOLD_PERCENTILE_SHORT_NAME);
        arguments.add("" + percentile);
        runCommandLine(arguments);
    }

    @Test(dataProvider="badPercentileMax50Data", expectedExceptions = UserException.BadArgumentValue.class)
    public void testBadOutlierColumnPercentile(final double percentile) {
        final List<String> arguments = new ArrayList<>();
        arguments.add("-" + StandardArgumentDefinitions.INPUT_SHORT_NAME);
        arguments.add(CONTROL_PCOV_FULL_FILE.toString());
        final File outputFile = createTempFile("pon-",".hd5");
        arguments.add("-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME);
        arguments.add(outputFile.toString());
        arguments.add("-" + CreatePanelOfNormals.COLUMN_EXTREME_THRESHOLD_PERCENTILE_SHORT_NAME);
        arguments.add("" + percentile);
        runCommandLine(arguments);
    }

    @Test(dataProvider="badNumOfEigenSamplesData", expectedExceptions = UserException.BadArgumentValue.class)
    public void testBadNumOfEigenSamples(final String numberOfEigenSamples) {
        final List<String> arguments = new ArrayList<>();
        arguments.add("-" + StandardArgumentDefinitions.INPUT_SHORT_NAME);
        arguments.add(CONTROL_PCOV_FULL_FILE.toString());
        final File outputFile = createTempFile("pon-",".hd5");
        arguments.add("-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME);
        arguments.add(outputFile.toString());
        arguments.add("-" + CreatePanelOfNormals.NUMBER_OF_EIGEN_SAMPLES_SHORT_NAME);
        arguments.add("" + numberOfEigenSamples);
        runCommandLine(arguments);
    }

    @DataProvider(name = "badNumOfEigenSamplesData")
    public Object[][] badNumOfEigenSamplesData() {
        return new Object[][]{
                { "-1.0" },
                { "" + Double.NaN },
                { "" + Double.POSITIVE_INFINITY },
                { "bad-number"},
                { "10.01" }
        };
    }


    @DataProvider(name = "badPercentileMax50Data")
    public Object[][] badPercentileMax50Data() {
        return new Object[][]{
                { -1.0 },
                { Double.NaN },
                { Double.POSITIVE_INFINITY },
                { 50.01 }
        };
    }

    @DataProvider(name = "badPercentileMax100Data")
    public Object[][] badPercentileMax100Data() {
        return new Object[][]{
                { -1.0 },
                { Double.NaN },
                { Double.POSITIVE_INFINITY },
                { 100.01 }
        };
    }

    @Test(dataProvider="someTargetsHDF5PoNCreationData")
    public void testSomeTargetsHDF5PoNCreation(final File targetsFile, final File inputFile) {
        final List<String> arguments = new ArrayList<>();
        arguments.add("-" + StandardArgumentDefinitions.INPUT_SHORT_NAME);
        arguments.add(inputFile.toString());
        if (targetsFile != null) {
            arguments.add("-" + TargetArgumentCollection.TARGET_FILE_SHORT_NAME);
            arguments.add(targetsFile.toString());
        }
        final File outputFile = createTempFile("pon-",".hd5");
        arguments.add("-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME);
        arguments.add(outputFile.toString());
        runCommandLine(arguments);
        assertEquivalentPoN(outputFile, EXPECTED_SOME_TARGETS_PON);
    }

    private void assertEquivalentPoN(final File outputFile, final File expectedAllTargetsPon) {
        try (final HDF5File leftFile = new HDF5File(outputFile);
             final HDF5File rightFile = new HDF5File(expectedAllTargetsPon)) {
            final HDF5PoN leftPoN = new HDF5PoN(leftFile);
            final HDF5PoN rightPoN = new HDF5PoN(rightFile);
            Assert.assertEquals(leftPoN.getSampleNames(), rightPoN.getSampleNames());
            Assert.assertEquals(new LinkedHashSet<>(leftPoN.getTargetNames()), new LinkedHashSet<>(rightPoN.getTargetNames()));
            assertEqualsMatrix(leftPoN.getTargetFactors(), rightPoN.getTargetFactors(), true);
            assertEqualsMatrix(leftPoN.getLogNormalizedCounts(), rightPoN.getLogNormalizedCounts(), true);
            assertEqualsMatrix(leftPoN.getLogNormalizedPInverseCounts(), rightPoN.getLogNormalizedPInverseCounts(), true);
            assertEqualsMatrix(leftPoN.getNormalizedCounts(), rightPoN.getNormalizedCounts(), true);
            assertEqualsMatrix(leftPoN.getReducedPanelCounts(), rightPoN.getReducedPanelCounts(), true);
            assertEqualsMatrix(leftPoN.getReducedPanelPInverseCounts(), rightPoN.getReducedPanelPInverseCounts(), true);
            Assert.assertEquals(leftPoN.getPanelSampleNames(), rightPoN.getPanelSampleNames());
            Assert.assertEquals(leftPoN.getPanelTargetNames(), rightPoN.getPanelTargetNames());
        }
    }

    private void assertEqualsMatrix(final RealMatrix left, final RealMatrix right, final boolean isAllowNegatedValues) {
        Assert.assertEquals(left.getRowDimension(), right.getRowDimension());
        Assert.assertEquals(left.getColumnDimension(), right.getColumnDimension());
        for (int i = 0; i < left.getRowDimension(); i++) {
            final double[] leftRow = left.getRow(i);
            final double[] rightRow = right.getRow(i);
            for (int j = 0; j < leftRow.length; j++) {
                if (isAllowNegatedValues) {
                    Assert.assertEquals(Math.abs(leftRow[j]), Math.abs(rightRow[j]), 0.0001);
                } else {
                    Assert.assertEquals(leftRow[j], rightRow[j], 0.0001);
                }
            }
        }
    }

    @DataProvider(name="allTargetsHDF5PoNCreationData")
    public Object[][] allTargetsHDF5PoNCreationData() {
        return new Object[][] {
                { null, CONTROL_PCOV_FULL_FILE },
                { ALL_TARGETS_FILE, CONTROL_PCOV_FULL_FILE },
                { ALL_TARGETS_FILE, CONTROL_PCOV_TARGET_NAME_ONLY_FILE },
                { ALL_TARGETS_FILE, CONTROL_PCOV_TARGET_COORDINATE_ONLY_FILE }
        };
    }

    @DataProvider(name="someTargetsHDF5PoNCreationData")
    public Object[][] someTargetsHDF5PoNCreationData() {
        return new Object[][] {
                { null, CONTROL_PCOV_SOME_TARGETS_FULL_FILE },
                { SOME_TARGETS_FILE, CONTROL_PCOV_FULL_FILE },
                { SOME_TARGETS_FILE, CONTROL_PCOV_TARGET_NAME_ONLY_FILE },
                { SOME_TARGETS_FILE, CONTROL_PCOV_TARGET_COORDINATE_ONLY_FILE },
                { SOME_TARGETS_FILE, CONTROL_PCOV_SOME_TARGETS_FULL_FILE },
        };
    }

}
