package org.broadinstitute.hellbender.tools.exome;

import org.apache.commons.io.FileUtils;
import org.apache.commons.lang.math.IntRange;
import org.broadinstitute.hdf5.HDF5File;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.hdf5.HDF5PoN;
import org.broadinstitute.hellbender.utils.hdf5.PoN;
import org.broadinstitute.hellbender.utils.hdf5.PoNTestUtils;
import org.broadinstitute.hellbender.utils.hdf5.RamPoN;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.IntStream;

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
    public void testAllTargetsHDF5PoNCreationNoSparkNoQC(final File targetsFile, final File inputFile) {
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
        arguments.add("--verbosity");
        arguments.add("INFO");
        arguments.add("-" + CreatePanelOfNormals.DISABLE_SPARK_SHORT_NAME);
        arguments.add("-" + CreatePanelOfNormals.NO_QC_SHORT_NAME);

        runCommandLine(arguments);
        assertBasicPoNAssumptions(outputFile, targetsFile);
        PoNTestUtils.assertEquivalentPoN(outputFile, EXPECTED_ALL_TARGETS_PON);
        assertRamPoNDuplicate(outputFile);

        // No blacklist file should have been created...
        final File blacklistFile = new File(outputFile + CreatePanelOfNormals.BLACKLIST_FILE_APPEND);
        Assert.assertFalse(blacklistFile.exists());
    }

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
        arguments.add("--verbosity");
        arguments.add("INFO");
        arguments.add("-" + CreatePanelOfNormals.DISABLE_SPARK_SHORT_NAME);

        runCommandLine(arguments);
        assertBasicPoNAssumptions(outputFile, targetsFile);
        PoNTestUtils.assertEquivalentPoN(outputFile, EXPECTED_ALL_TARGETS_PON);
        assertRamPoNDuplicate(outputFile);

        final File blacklistFile = new File(outputFile + CreatePanelOfNormals.BLACKLIST_FILE_APPEND);
        assertBlacklistFileIsEmpty(blacklistFile);
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
        final File outputFile = createTempFile("pon-", ".hd5");
        arguments.add("-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME);
        arguments.add(outputFile.toString());
        arguments.add("--verbosity");
        arguments.add("INFO");

        runCommandLine(arguments);
        // Need to support normalizedCoveragedFile == null
        assertBasicPoNAssumptions(outputFile, targetsFile);
        PoNTestUtils.assertEquivalentPoN(outputFile, EXPECTED_ALL_TARGETS_PON);
        assertRamPoNDuplicate(outputFile);

        final File blacklistFile = new File(outputFile + CreatePanelOfNormals.BLACKLIST_FILE_APPEND);
        assertBlacklistFileIsEmpty(blacklistFile);
    }

    @Test(dataProvider="allTargetsHDF5PoNCreationData")
    public void testAllTargetsHDF5PoNCreationSparkNoQC(final File targetsFile, final File inputFile) {
        final List<String> arguments = new ArrayList<>();
        arguments.add("-" + StandardArgumentDefinitions.INPUT_SHORT_NAME);
        arguments.add(inputFile.toString());
        if (targetsFile != null) {
            arguments.add("-" + TargetArgumentCollection.TARGET_FILE_SHORT_NAME);
            arguments.add(targetsFile.toString());
        }
        final File outputFile = createTempFile("pon-", ".hd5");
        arguments.add("-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME);
        arguments.add(outputFile.toString());
        arguments.add("--verbosity");
        arguments.add("INFO");
        arguments.add("-" + CreatePanelOfNormals.NO_QC_SHORT_NAME);

        runCommandLine(arguments);
        // Need to support normalizedCoveragedFile == null
        assertBasicPoNAssumptions(outputFile, targetsFile);
        PoNTestUtils.assertEquivalentPoN(outputFile, EXPECTED_ALL_TARGETS_PON);
        assertRamPoNDuplicate(outputFile);

        // No blacklist file should have been created...
        final File blacklistFile = new File(outputFile + CreatePanelOfNormals.BLACKLIST_FILE_APPEND);
        Assert.assertFalse(blacklistFile.exists());
    }

    private void assertBlacklistFileIsEmpty(final File blacklistFile) {
        Assert.assertTrue(blacklistFile.exists());
        Assert.assertEquals(FileUtils.sizeOf(blacklistFile), 0);
    }

    private void assertBasicPoNAssumptions(final File ponFile, final File initialTargetsFileUsedToCreatePoN) {
        try (final HDF5File ponHDF5File = new HDF5File(ponFile)) {
            final HDF5PoN pon = new HDF5PoN(ponHDF5File);

            Assert.assertTrue(pon.getTargets().size() >= pon.getPanelTargets().size());
            Assert.assertTrue(pon.getRawTargets().size() > pon.getTargets().size());

            Assert.assertTrue(pon.getTargetNames().size() == pon.getTargets().size());
            Assert.assertTrue(pon.getPanelTargetNames().size() == pon.getPanelTargetNames().size());
            Assert.assertTrue(pon.getRawTargetNames().size() == pon.getRawTargetNames().size());

            if (initialTargetsFileUsedToCreatePoN != null) {

                final TargetCollection<Target> tc = TargetArgumentCollection.readTargetCollection(initialTargetsFileUsedToCreatePoN);
                Assert.assertEquals(pon.getRawTargets().size(), tc.targetCount());

                // Check that the raw targets are the same
                Assert.assertTrue(IntStream.of(new IntRange(0, pon.getRawTargets().size()-1).toArray()).boxed().map(i -> pon.getRawTargets().get(i).equals(tc.target(i))).allMatch(t -> t));
            }
        }
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
    public void testSomeTargetsHDF5PoNCreationNoSpark(final File targetsFile, final File inputFile) {
        final List<String> arguments = new ArrayList<>();
        arguments.add("-" + StandardArgumentDefinitions.INPUT_SHORT_NAME);
        arguments.add(inputFile.toString());
        if (targetsFile != null) {
            arguments.add("-" + TargetArgumentCollection.TARGET_FILE_SHORT_NAME);
            arguments.add(targetsFile.toString());
        }
        arguments.add("-ds");
        final File outputFile = createTempFile("pon-",".hd5");
        arguments.add("-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME);
        arguments.add(outputFile.toString());
        runCommandLine(arguments);
        PoNTestUtils.assertEquivalentPoN(outputFile, EXPECTED_SOME_TARGETS_PON);
        assertRamPoNDuplicate(outputFile);
    }

    @Test(dataProvider="someTargetsHDF5PoNCreationData")
    public void testSomeTargetsHDF5PoNCreationSpark(final File targetsFile, final File inputFile) {
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
        PoNTestUtils.assertEquivalentPoN(outputFile, EXPECTED_SOME_TARGETS_PON);
        assertRamPoNDuplicate(outputFile);
    }


    private void assertRamPoNDuplicate(final File outputFile) {
        try (final HDF5File hdf5FilePoN = new HDF5File(outputFile)) {
            final HDF5PoN filePoN = new HDF5PoN(hdf5FilePoN);
            assertRamPoNDuplicate(filePoN);
        }
    }

    /**
     * Assert that we can create a valid ram pon duplicate.
     *
     * @param pon never {@code null}
     */
    private void assertRamPoNDuplicate(final PoN pon) {

        final PoN ramPoN = new RamPoN(pon);

        PoNTestUtils.assertEqualsMatrix(pon.getLogNormalizedCounts(), ramPoN.getLogNormalizedCounts(), false);
        PoNTestUtils.assertEqualsMatrix(pon.getLogNormalizedPInverseCounts(), ramPoN.getLogNormalizedPInverseCounts(), false);
        PoNTestUtils.assertEqualsMatrix(pon.getNormalizedCounts(), ramPoN.getNormalizedCounts(), false);
        PoNTestUtils.assertEqualsMatrix(pon.getReducedPanelCounts(), ramPoN.getReducedPanelCounts(), false);
        PoNTestUtils.assertEqualsMatrix(pon.getReducedPanelPInverseCounts(), ramPoN.getReducedPanelPInverseCounts(), false);
        PoNTestUtils.assertEqualsMatrix(pon.getTargetFactors(), ramPoN.getTargetFactors(), false);

        final List<Target> ponPanelTargets = pon.getPanelTargets();
        final List<Target> ponTargets = pon.getTargets();
        final List<Target> ponRawTargets = pon.getRawTargets();

        final List<Target> ramPoNPanelTargets = ramPoN.getPanelTargets();
        final List<Target> ramPoNTargets = ramPoN.getTargets();
        final List<Target> ramPoNRawTargets = ramPoN.getRawTargets();

        Assert.assertEquals(ponPanelTargets.size(), ramPoNPanelTargets.size());
        Assert.assertEquals(ponTargets.size(), ramPoNTargets.size());

        // Make sure every target is the same
        Assert.assertTrue(IntStream.range(0, ponRawTargets.size()).boxed().allMatch(i -> ponRawTargets.get(i).equals(ramPoNRawTargets.get(i))));
        Assert.assertTrue(IntStream.range(0, ponTargets.size()).boxed().allMatch(i -> ponTargets.get(i).equals(ramPoNTargets.get(i))));
        Assert.assertTrue(IntStream.range(0, ponPanelTargets.size()).boxed().allMatch(i -> ponPanelTargets.get(i).equals(ramPoNPanelTargets.get(i))));

        // Make sure every sample name is the same
        Assert.assertTrue(IntStream.range(0, pon.getSampleNames().size()).boxed().allMatch(i -> pon.getSampleNames().get(i).equals(ramPoN.getSampleNames().get(i))));
        Assert.assertTrue(IntStream.range(0, pon.getPanelSampleNames().size()).boxed().allMatch(i -> pon.getPanelSampleNames().get(i).equals(ramPoN.getPanelSampleNames().get(i))));

        // Make sure every target name is the same
        Assert.assertTrue(IntStream.range(0, pon.getTargetNames().size()).boxed().allMatch(i -> pon.getTargetNames().get(i).equals(ramPoN.getTargetNames().get(i))));
        Assert.assertTrue(IntStream.range(0, pon.getPanelTargetNames().size()).boxed().allMatch(i -> pon.getPanelTargetNames().get(i).equals(ramPoN.getPanelTargetNames().get(i))));

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
