package org.broadinstitute.hellbender.tools.exome;

import htsjdk.tribble.bed.BEDCodec;
import htsjdk.tribble.bed.BEDFeature;
import org.apache.commons.math3.linear.RealMatrix;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.hdf5.HDF5LibraryUnitTests;
import org.broadinstitute.hellbender.utils.hdf5.HDF5PoN;
import org.broadinstitute.hellbender.utils.hdf5.HDF5Reader;
import org.broadinstitute.hellbender.utils.hdf5.PoN;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

/**
 * Integration tests for {@link NormalizeSomaticReadCounts}.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public class NormalizeSomaticReadCountsIntegrationTest extends CommandLineProgramTest {


    private static final File TEST_DIR = new File("src/test/resources/org/broadinstitute/tools/exome");

    private static final File FULL_READ_COUNTS_INPUT = new File(TEST_DIR,"full-read-counts.txt");
    private static final File TARGET_NAME_ONLY_READ_COUNTS_INPUT = new File(TEST_DIR,"only-names-read-counts.txt");
    private static final File COORD_ONLY_READ_COUNTS_INPUT = new File(TEST_DIR,"only-coords-read-counts.txt");
    private static final File FULL_READ_COUNTS_WITH_EXTRA_TARGET_INPUT = new File(TEST_DIR,"full-read-counts-with-extra-target.txt");
    private static final File FULL_READ_COUNTS_MISSING_A_TARGET_INPUT = new File(TEST_DIR,"full-read-counts-missing-a-target.txt");

    private static final File TEST_TARGETS = new File(TEST_DIR,"targets.bed");

    private static final File TEST_TARGETS_WITH_BAD_NAME = new File(TEST_DIR,"targets-with-bad-name.bed");

    private static final File TEST_PON = HDF5LibraryUnitTests.TEST_PON;


    @Override
    public String getTestedClassName() {
        return NormalizeSomaticReadCounts.class.getSimpleName();
    }

    @Test(expectedExceptions = UserException.BadInput.class)
    public void testRunWithTargetFileWithBadName() throws IOException {
        final File factorNormalizedOutput = createTempFile("test",".txt");

        final String[] arguments = {
                "-" + NormalizeSomaticReadCounts.READ_COUNTS_FILE_SHORT_NAME, COORD_ONLY_READ_COUNTS_INPUT.getAbsolutePath(),
                "-" + NormalizeSomaticReadCounts.PON_FILE_SHORT_NAME, TEST_PON.getAbsolutePath(),
                "-" + NormalizeSomaticReadCounts.FACTOR_NORMALIZED_COUNTS_SHORT_NAME, factorNormalizedOutput.getAbsolutePath(),
                "-" + NormalizeSomaticReadCounts.TARGET_FILE_SHORT_NAME, TEST_TARGETS_WITH_BAD_NAME.getAbsolutePath(),
        };
        runCommandLine(arguments);
    }

    @Test(expectedExceptions = UserException.CouldNotReadInputFile.class)
    public void testBadTargetFile() throws IOException {
        final File factorNormalizedOutput = createTempFile("test",".txt");

        final String[] arguments = {
                "-" + NormalizeSomaticReadCounts.READ_COUNTS_FILE_SHORT_NAME, FULL_READ_COUNTS_INPUT.getAbsolutePath(),
                "-" + NormalizeSomaticReadCounts.PON_FILE_SHORT_NAME, TEST_PON.getAbsolutePath(),
                "-" + NormalizeSomaticReadCounts.FACTOR_NORMALIZED_COUNTS_SHORT_NAME, factorNormalizedOutput.getAbsolutePath(),
                "-" + NormalizeSomaticReadCounts.TARGET_FILE_SHORT_NAME, TEST_TARGETS_WITH_BAD_NAME.getAbsolutePath() + "failure-name",
        };
        runCommandLine(arguments);
    }

    @Test
    public void testFullInputWithExtraTarget() throws IOException {
        final File factorNormalizedOutput = createTempFile("test",".txt");

        final String[] arguments = {
                "-" + NormalizeSomaticReadCounts.READ_COUNTS_FILE_SHORT_NAME, FULL_READ_COUNTS_WITH_EXTRA_TARGET_INPUT.getAbsolutePath(),
                "-" + NormalizeSomaticReadCounts.PON_FILE_SHORT_NAME, TEST_PON.getAbsolutePath(),
                "-" + NormalizeSomaticReadCounts.FACTOR_NORMALIZED_COUNTS_SHORT_NAME, factorNormalizedOutput.getAbsolutePath(),
        };
        runCommandLine(arguments);
    }

    @Test(expectedExceptions = UserException.BadInput.class)
    public void testFullInputMissingATarget() throws IOException {
        final File factorNormalizedOutput = createTempFile("test",".txt");

        final String[] arguments = {
                "-" + NormalizeSomaticReadCounts.READ_COUNTS_FILE_SHORT_NAME, FULL_READ_COUNTS_MISSING_A_TARGET_INPUT.getAbsolutePath(),
                "-" + NormalizeSomaticReadCounts.PON_FILE_SHORT_NAME, TEST_PON.getAbsolutePath(),
                "-" + NormalizeSomaticReadCounts.FACTOR_NORMALIZED_COUNTS_SHORT_NAME, factorNormalizedOutput.getAbsolutePath(),
        };
        runCommandLine(arguments);
    }


    @Test
    public void testFullInputRunWithTargetFileWithBadName() throws IOException {
        final File factorNormalizedOutput = createTempFile("test",".txt");

        final String[] arguments = {
                "-" + NormalizeSomaticReadCounts.READ_COUNTS_FILE_SHORT_NAME, FULL_READ_COUNTS_INPUT.getAbsolutePath(),
                "-" + NormalizeSomaticReadCounts.PON_FILE_SHORT_NAME, TEST_PON.getAbsolutePath(),
                "-" + NormalizeSomaticReadCounts.FACTOR_NORMALIZED_COUNTS_SHORT_NAME, factorNormalizedOutput.getAbsolutePath(),
                "-" + NormalizeSomaticReadCounts.TARGET_FILE_SHORT_NAME, TEST_TARGETS_WITH_BAD_NAME.getAbsolutePath(),
        };
        runCommandLine(arguments);
    }

    @Test
    public void testFullReadCountsInputRun() throws IOException {
       final File factorNormalizedOutput = createTempFile("test",".txt");

       final String[] arguments = {
               "-" + NormalizeSomaticReadCounts.READ_COUNTS_FILE_SHORT_NAME, FULL_READ_COUNTS_INPUT.getAbsolutePath(),
               "-" + NormalizeSomaticReadCounts.PON_FILE_SHORT_NAME, TEST_PON.getAbsolutePath(),
               "-" + NormalizeSomaticReadCounts.FACTOR_NORMALIZED_COUNTS_SHORT_NAME, factorNormalizedOutput.getAbsolutePath()
       };

       runCommandLine(arguments);
       final ReadCountCollection input = ReadCountCollectionUtils.parse(FULL_READ_COUNTS_INPUT, null);
       final ReadCountCollection factorNormalized = ReadCountCollectionUtils.parse(factorNormalizedOutput, null);

       Assert.assertEquals(factorNormalized.columnNames(), input.columnNames());
       Assert.assertTrue(!factorNormalized.targets().stream().anyMatch(t -> t.getInterval() == null));
       Assert.assertEquals(factorNormalized.targets().stream().map(Target::getInterval).collect(Collectors.toSet()),
               input.targets().stream().map(Target::getInterval).collect(Collectors.toSet()));
       Assert.assertEquals(factorNormalized.targets().stream().collect(Collectors.toSet()),
               input.targets().stream().collect(Collectors.toSet()));
       assertFactorNormalizedValues(input, factorNormalized);
    }

    @Test
    public void testNameOnlyCountsInputRun() throws IOException {
        final File factorNormalizedOutput = createTempFile("test",".txt");

        final String[] arguments = {
                "-" + NormalizeSomaticReadCounts.READ_COUNTS_FILE_SHORT_NAME, TARGET_NAME_ONLY_READ_COUNTS_INPUT.getAbsolutePath(),
                "-" + NormalizeSomaticReadCounts.PON_FILE_SHORT_NAME, TEST_PON.getAbsolutePath(),
                "-" + NormalizeSomaticReadCounts.FACTOR_NORMALIZED_COUNTS_SHORT_NAME, factorNormalizedOutput.getAbsolutePath(),
        };

        runCommandLine(arguments);
        final ReadCountCollection input = ReadCountCollectionUtils.parse(TARGET_NAME_ONLY_READ_COUNTS_INPUT, null);
        final ReadCountCollection factorNormalized = ReadCountCollectionUtils.parse(factorNormalizedOutput, null);

        Assert.assertEquals(factorNormalized.columnNames(), input.columnNames());
        Assert.assertFalse(factorNormalized.targets().stream().anyMatch(t -> t.getInterval() != null));
        Assert.assertEquals(factorNormalized.targets().stream().collect(Collectors.toSet()),
                input.targets().stream().collect(Collectors.toSet()));
        assertFactorNormalizedValues(input,factorNormalized);
    }

    @Test
    public void testCoordOnlyCountsInputFileRun() throws IOException {
        final File factorNormalizedOutput = createTempFile("test",".txt");

        final String[] arguments = {
                "-" + NormalizeSomaticReadCounts.READ_COUNTS_FILE_SHORT_NAME, COORD_ONLY_READ_COUNTS_INPUT.getAbsolutePath(),
                "-" + NormalizeSomaticReadCounts.PON_FILE_SHORT_NAME, TEST_PON.getAbsolutePath(),
                "-" + NormalizeSomaticReadCounts.FACTOR_NORMALIZED_COUNTS_SHORT_NAME, factorNormalizedOutput.getAbsolutePath(),
                "-" + NormalizeSomaticReadCounts.TARGET_FILE_SHORT_NAME, TEST_TARGETS.getAbsolutePath(),
        };

        runCommandLine(arguments);
        final ExonCollection<? extends BEDFeature> exons = ExonCollections.fromBEDFeatureFile(TEST_TARGETS,new BEDCodec());
        final ReadCountCollection input = ReadCountCollectionUtils.parse(COORD_ONLY_READ_COUNTS_INPUT,exons);
        final ReadCountCollection factorNormalized = ReadCountCollectionUtils.parse(factorNormalizedOutput, exons);

        Assert.assertEquals(factorNormalized.columnNames(), input.columnNames());
        Assert.assertFalse(factorNormalized.targets().stream().anyMatch(t -> t.getInterval() == null));
       // Assert.assertTrue(factorNormalized.hasTargetNames());
       // Assert.assertEquals(factorNormalized.getIntervals(), input.getIntervals());
        assertFactorNormalizedValues(input, factorNormalized);
    }

    @Test(expectedExceptions = UserException.BadInput.class)
    public void testCoordOnlyCountsMissingTargetInputFileRun() throws IOException {
        final File factorNormalizedOutput = createTempFile("test",".txt");

        final String[] arguments = {
                "-" + NormalizeSomaticReadCounts.READ_COUNTS_FILE_SHORT_NAME, COORD_ONLY_READ_COUNTS_INPUT.getAbsolutePath(),
                "-" + NormalizeSomaticReadCounts.PON_FILE_SHORT_NAME, TEST_PON.getAbsolutePath(),
                "-" + NormalizeSomaticReadCounts.FACTOR_NORMALIZED_COUNTS_SHORT_NAME, factorNormalizedOutput.getAbsolutePath(),
        };

        runCommandLine(arguments);
    }

    @DataProvider(name="inputFileData")
    public Object[][] inputFileData() {
        return new Object[][] {
                new Object[] { FULL_READ_COUNTS_INPUT },
                new Object[] { TARGET_NAME_ONLY_READ_COUNTS_INPUT },
                new Object[] { COORD_ONLY_READ_COUNTS_INPUT },
        };
    }

    private void assertFactorNormalizedValues(final ReadCountCollection input, final ReadCountCollection factorNormalized) {
        try (final HDF5Reader ponReader = new HDF5Reader(TEST_PON)) {
            final PoN pon = new HDF5PoN(ponReader);
            final RealMatrix targetFactors = pon.targetFactors();
            final List<String> ponTargets = pon.targetNames();
            final Map<String,Integer> ponTargetIndexes = new HashMap<>(ponTargets.size());
            for (int i = 0; i < ponTargets.size(); i++) {
                ponTargetIndexes.put(ponTargets.get(i),i);
            }
            final RealMatrix inputCounts = input.counts();
            final RealMatrix factorNormalizedCounts = factorNormalized.counts();
            for (int i = 0; i < factorNormalizedCounts.getRowDimension(); i++) {
                final double factor = targetFactors.getEntry(ponTargetIndexes.get(factorNormalized.targets().get(i).getName()),0);
                final double[] inputValues = inputCounts.getRow(i);
                final double[] outputValues = factorNormalizedCounts.getRow(i);
                for (int j = 0; j < inputValues.length; j++) {
                    final double expected = inputValues[j] / factor;
                    Assert.assertEquals(outputValues[j],expected,0.0000001,"" + i + " , " + j);
                }
            }
        }
    }
}
