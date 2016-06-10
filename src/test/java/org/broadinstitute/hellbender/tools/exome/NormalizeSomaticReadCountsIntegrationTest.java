package org.broadinstitute.hellbender.tools.exome;

import htsjdk.tribble.bed.BEDCodec;
import htsjdk.tribble.bed.BEDFeature;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.DefaultRealMatrixChangingVisitor;
import org.apache.commons.math3.linear.DefaultRealMatrixPreservingVisitor;
import org.apache.commons.math3.linear.RealMatrix;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.ExomeStandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.GATKProtectedMathUtils;
import org.broadinstitute.hellbender.utils.hdf5.HDF5File;
import org.broadinstitute.hellbender.utils.hdf5.HDF5LibraryUnitTest;
import org.broadinstitute.hellbender.utils.hdf5.HDF5PoN;
import org.broadinstitute.hellbender.utils.hdf5.PoN;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableReader;
import org.broadinstitute.hellbender.utils.tsv.TableUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;

/**
 * Integration tests for {@link NormalizeSomaticReadCounts}.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public class NormalizeSomaticReadCountsIntegrationTest extends CommandLineProgramTest {

    private static final File TEST_DIR = new File("src/test/resources/org/broadinstitute/hellbender/tools/exome");

    private static final File FULL_READ_COUNTS_INPUT = new File(TEST_DIR,"full-read-counts.txt");
    private static final File FULL_READ_COUNTS_INPUT_ONE_SAMPLE = new File(TEST_DIR,"full-read-counts.1sample.txt");
    private static final File TARGET_NAME_ONLY_READ_COUNTS_INPUT = new File(TEST_DIR,"only-names-read-counts.txt");
    private static final File TARGET_NAME_ONLY_READ_COUNTS_INPUT_ONE_SAMPLE = new File(TEST_DIR,"only-names-read-counts.1sample.txt");
    private static final File COORD_ONLY_READ_COUNTS_INPUT = new File(TEST_DIR,"only-coords-read-counts.txt");
    private static final File COORD_ONLY_READ_COUNTS_INPUT_ONE_SAMPLE = new File(TEST_DIR,"only-coords-read-counts.1sample.txt");
    private static final File FULL_READ_COUNTS_WITH_EXTRA_TARGET_INPUT = new File(TEST_DIR,"full-read-counts-with-extra-target.txt");
    private static final File FULL_READ_COUNTS_WITH_EXTRA_TARGET_INPUT_ONE_SAMPLE = new File(TEST_DIR,"full-read-counts-with-extra-target.1sample.txt");
    private static final File FULL_READ_COUNTS_MISSING_A_TARGET_INPUT = new File(TEST_DIR,"full-read-counts-missing-a-target.txt");
    private static final File FULL_READ_COUNTS_MISSING_A_TARGET_INPUT_ONE_SAMPLE = new File(TEST_DIR,"full-read-counts-missing-a-target.1sample.txt");
    private static final File FULL_READ_COUNTS_BAD_NAME = new File(TEST_DIR,"full-read-counts-bad-target-name.txt");
    private static final File FULL_READ_COUNTS_BAD_NAME_ONE_SAMPLE = new File(TEST_DIR,"full-read-counts-bad-target-name.1sample.txt");

    private static final File TEST_TARGETS = new File(TEST_DIR,"targets.bed");

    private static final File TEST_TARGETS_WITH_BAD_NAME = new File(TEST_DIR,"targets-with-bad-name.bed");

    private static final File TEST_PON = HDF5LibraryUnitTest.TEST_PON;


    @Override
    public String getTestedClassName() {
        return NormalizeSomaticReadCounts.class.getSimpleName();
    }

    @Test(expectedExceptions = UserException.BadInput.class)
    public void testRunWithTargetFileWithBadName() throws IOException {
        final File factorNormalizedOutput = createTempFile("test",".txt");
        final File tangentNormalizationOutput = createTempFile("test",".txt");

        final String[] arguments = {
                "-" + NormalizeSomaticReadCounts.READ_COUNTS_FILE_SHORT_NAME, COORD_ONLY_READ_COUNTS_INPUT.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.PON_FILE_SHORT_NAME, TEST_PON.getAbsolutePath(),
                "-" + NormalizeSomaticReadCounts.FACTOR_NORMALIZED_COUNTS_SHORT_NAME, factorNormalizedOutput.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.TANGENT_NORMALIZED_COUNTS_FILE_SHORT_NAME, tangentNormalizationOutput.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.TARGET_FILE_SHORT_NAME, TEST_TARGETS_WITH_BAD_NAME.getAbsolutePath(),
        };
        runCommandLine(arguments);
    }

    @Test(expectedExceptions = UserException.CouldNotReadInputFile.class)
    public void testBadTargetFile() throws IOException {
        final File factorNormalizedOutput = createTempFile("test",".txt");
        final File tangentNormalizationOutput = createTempFile("test",".txt");
        final File betaHatsOutput = createTempFile("tangent-", ".bhats");
        final File preTangentNormalizationOutput = createTempFile("pre-tn-",".txt");

        final String[] arguments = {
                "-" + NormalizeSomaticReadCounts.READ_COUNTS_FILE_SHORT_NAME, FULL_READ_COUNTS_INPUT.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.PON_FILE_SHORT_NAME, TEST_PON.getAbsolutePath(),
                "-" + NormalizeSomaticReadCounts.FACTOR_NORMALIZED_COUNTS_SHORT_NAME, factorNormalizedOutput.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.TANGENT_NORMALIZED_COUNTS_FILE_SHORT_NAME, tangentNormalizationOutput.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.TARGET_FILE_SHORT_NAME, TEST_TARGETS_WITH_BAD_NAME.getAbsolutePath() + "failure-name",
                "-" + NormalizeSomaticReadCounts.TANGENT_BETA_HATS_SHORT_NAME, betaHatsOutput.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.PRE_TANGENT_NORMALIZED_COUNTS_FILE_SHORT_NAME, preTangentNormalizationOutput.getAbsolutePath()
        };
        runCommandLine(arguments);
    }

    @Test(expectedExceptions = UserException.BadInput.class)
    public void testFullInputWithExtraTargetMultipleSamples() throws IOException {
        final File factorNormalizedOutput = createTempFile("test",".txt");
        final File tangentNormalizationOutput = createTempFile("test",".txt");
        final File betaHatsOutput = createTempFile("tangent-", ".bhats");
        final File preTangentNormalizationOutput = createTempFile("pre-tn-",".txt");

        final String[] arguments = {
                "-" + NormalizeSomaticReadCounts.READ_COUNTS_FILE_SHORT_NAME, FULL_READ_COUNTS_WITH_EXTRA_TARGET_INPUT.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.PON_FILE_SHORT_NAME, TEST_PON.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.TANGENT_NORMALIZED_COUNTS_FILE_SHORT_NAME, tangentNormalizationOutput.getAbsolutePath(),
                "-" + NormalizeSomaticReadCounts.FACTOR_NORMALIZED_COUNTS_SHORT_NAME, factorNormalizedOutput.getAbsolutePath(),
                "-" + NormalizeSomaticReadCounts.TANGENT_BETA_HATS_SHORT_NAME, betaHatsOutput.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.PRE_TANGENT_NORMALIZED_COUNTS_FILE_SHORT_NAME, preTangentNormalizationOutput.getAbsolutePath()
        };
        runCommandLine(arguments);
    }

    @Test
    public void testFullInputWithExtraTargetOneSample() throws IOException {
        final File factorNormalizedOutput = createTempFile("test",".txt");
        final File tangentNormalizationOutput = createTempFile("test",".txt");
        final File betaHatsOutput = createTempFile("tangent-", ".bhats");
        final File preTangentNormalizationOutput = createTempFile("pre-tn-",".txt");

        final String[] arguments = {
                "-" + NormalizeSomaticReadCounts.READ_COUNTS_FILE_SHORT_NAME, FULL_READ_COUNTS_WITH_EXTRA_TARGET_INPUT_ONE_SAMPLE.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.PON_FILE_SHORT_NAME, TEST_PON.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.TANGENT_NORMALIZED_COUNTS_FILE_SHORT_NAME, tangentNormalizationOutput.getAbsolutePath(),
                "-" + NormalizeSomaticReadCounts.FACTOR_NORMALIZED_COUNTS_SHORT_NAME, factorNormalizedOutput.getAbsolutePath(),
                "-" + NormalizeSomaticReadCounts.TANGENT_BETA_HATS_SHORT_NAME, betaHatsOutput.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.PRE_TANGENT_NORMALIZED_COUNTS_FILE_SHORT_NAME, preTangentNormalizationOutput.getAbsolutePath()
        };
        runCommandLine(arguments);
    }

    @Test(expectedExceptions = UserException.BadInput.class)
    public void testFullInputMissingATargetMultipleSamples() throws IOException {
        final File factorNormalizedOutput = createTempFile("test",".txt");
        final File tangentNormalizationOutput = createTempFile("test",".txt");
        final File betaHatsOutput = createTempFile("tangent-", ".bhats");
        final File preTangentNormalizationOutput = createTempFile("pre-tn-",".txt");


        final String[] arguments = {
                "-" + NormalizeSomaticReadCounts.READ_COUNTS_FILE_SHORT_NAME, FULL_READ_COUNTS_MISSING_A_TARGET_INPUT.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.PON_FILE_SHORT_NAME, TEST_PON.getAbsolutePath(),
                "-" + NormalizeSomaticReadCounts.FACTOR_NORMALIZED_COUNTS_SHORT_NAME, factorNormalizedOutput.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.TANGENT_NORMALIZED_COUNTS_FILE_SHORT_NAME, tangentNormalizationOutput.getAbsolutePath(),
                "-" + NormalizeSomaticReadCounts.TANGENT_BETA_HATS_SHORT_NAME, betaHatsOutput.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.PRE_TANGENT_NORMALIZED_COUNTS_FILE_SHORT_NAME, preTangentNormalizationOutput.getAbsolutePath()
        };
        runCommandLine(arguments);
    }

    @Test(expectedExceptions = UserException.BadInput.class)
    public void testFullInputMissingATargetOneSample() throws IOException {
        final File factorNormalizedOutput = createTempFile("test",".txt");
        final File tangentNormalizationOutput = createTempFile("test",".txt");
        final File betaHatsOutput = createTempFile("tangent-", ".bhats");
        final File preTangentNormalizationOutput = createTempFile("pre-tn-",".txt");


        final String[] arguments = {
                "-" + NormalizeSomaticReadCounts.READ_COUNTS_FILE_SHORT_NAME, FULL_READ_COUNTS_MISSING_A_TARGET_INPUT_ONE_SAMPLE.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.PON_FILE_SHORT_NAME, TEST_PON.getAbsolutePath(),
                "-" + NormalizeSomaticReadCounts.FACTOR_NORMALIZED_COUNTS_SHORT_NAME, factorNormalizedOutput.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.TANGENT_NORMALIZED_COUNTS_FILE_SHORT_NAME, tangentNormalizationOutput.getAbsolutePath(),
                "-" + NormalizeSomaticReadCounts.TANGENT_BETA_HATS_SHORT_NAME, betaHatsOutput.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.PRE_TANGENT_NORMALIZED_COUNTS_FILE_SHORT_NAME, preTangentNormalizationOutput.getAbsolutePath()
        };
        runCommandLine(arguments);
    }


    @Test(expectedExceptions = UserException.BadInput.class)
    public void testFullInputRunWithTargetFileWithBadNameMultipleSamples() throws IOException {
        final File factorNormalizedOutput = createTempFile("test",".txt");
        final File tangentNormalizationOutput = createTempFile("test",".txt");


        final String[] arguments = {
                "-" + NormalizeSomaticReadCounts.READ_COUNTS_FILE_SHORT_NAME, FULL_READ_COUNTS_BAD_NAME.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.PON_FILE_SHORT_NAME, TEST_PON.getAbsolutePath(),
                "-" + NormalizeSomaticReadCounts.FACTOR_NORMALIZED_COUNTS_SHORT_NAME, factorNormalizedOutput.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.TARGET_FILE_SHORT_NAME, TEST_TARGETS_WITH_BAD_NAME.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.TANGENT_NORMALIZED_COUNTS_FILE_SHORT_NAME, tangentNormalizationOutput.getAbsolutePath(),
        };
        runCommandLine(arguments);
    }

    @Test(expectedExceptions = UserException.BadInput.class)
    public void testFullInputRunWithTargetFileWithBadNameOneSample() throws IOException {
        final File factorNormalizedOutput = createTempFile("test",".txt");
        final File tangentNormalizationOutput = createTempFile("test",".txt");


        final String[] arguments = {
                "-" + NormalizeSomaticReadCounts.READ_COUNTS_FILE_SHORT_NAME, FULL_READ_COUNTS_BAD_NAME_ONE_SAMPLE.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.PON_FILE_SHORT_NAME, TEST_PON.getAbsolutePath(),
                "-" + NormalizeSomaticReadCounts.FACTOR_NORMALIZED_COUNTS_SHORT_NAME, factorNormalizedOutput.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.TARGET_FILE_SHORT_NAME, TEST_TARGETS_WITH_BAD_NAME.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.TANGENT_NORMALIZED_COUNTS_FILE_SHORT_NAME, tangentNormalizationOutput.getAbsolutePath(),
        };
        runCommandLine(arguments);
    }

    @Test(expectedExceptions = UserException.BadInput.class)
    public void testFullReadCountsInputRunMultipleSamples() throws IOException {
        final File factorNormalizedOutput = createTempFile("test",".txt");
        final File tangentNormalizationOutput = createTempFile("test",".txt");
        final File betaHatsOutput = createTempFile("tangent-", ".bhats");
        final File preTangentNormalizationOutput = createTempFile("pre-tn-",".txt");


        final String[] arguments = {
               "-" + NormalizeSomaticReadCounts.READ_COUNTS_FILE_SHORT_NAME, FULL_READ_COUNTS_INPUT.getAbsolutePath(),
               "-" + ExomeStandardArgumentDefinitions.PON_FILE_SHORT_NAME, TEST_PON.getAbsolutePath(),
               "-" + NormalizeSomaticReadCounts.FACTOR_NORMALIZED_COUNTS_SHORT_NAME, factorNormalizedOutput.getAbsolutePath(),
               "-" + ExomeStandardArgumentDefinitions.TANGENT_NORMALIZED_COUNTS_FILE_SHORT_NAME, tangentNormalizationOutput.getAbsolutePath(),
               "-" + NormalizeSomaticReadCounts.TANGENT_BETA_HATS_SHORT_NAME, betaHatsOutput.getAbsolutePath(),
               "-" + ExomeStandardArgumentDefinitions.PRE_TANGENT_NORMALIZED_COUNTS_FILE_SHORT_NAME, preTangentNormalizationOutput.getAbsolutePath()
        };

        runCommandLine(arguments);
    }

    @Test
    public void testFullReadCountsInputRunOneSample() throws IOException {
        final File factorNormalizedOutput = createTempFile("test",".txt");
        final File tangentNormalizationOutput = createTempFile("test",".txt");
        final File betaHatsOutput = createTempFile("tangent-", ".bhats");
        final File preTangentNormalizationOutput = createTempFile("pre-tn-",".txt");


        final String[] arguments = {
               "-" + NormalizeSomaticReadCounts.READ_COUNTS_FILE_SHORT_NAME, FULL_READ_COUNTS_INPUT_ONE_SAMPLE.getAbsolutePath(),
               "-" + ExomeStandardArgumentDefinitions.PON_FILE_SHORT_NAME, TEST_PON.getAbsolutePath(),
               "-" + NormalizeSomaticReadCounts.FACTOR_NORMALIZED_COUNTS_SHORT_NAME, factorNormalizedOutput.getAbsolutePath(),
               "-" + ExomeStandardArgumentDefinitions.TANGENT_NORMALIZED_COUNTS_FILE_SHORT_NAME, tangentNormalizationOutput.getAbsolutePath(),
               "-" + NormalizeSomaticReadCounts.TANGENT_BETA_HATS_SHORT_NAME, betaHatsOutput.getAbsolutePath(),
               "-" + ExomeStandardArgumentDefinitions.PRE_TANGENT_NORMALIZED_COUNTS_FILE_SHORT_NAME, preTangentNormalizationOutput.getAbsolutePath()
        };

        runCommandLine(arguments);
        final ReadCountCollection input = ReadCountCollectionUtils.parse(FULL_READ_COUNTS_INPUT_ONE_SAMPLE);
        final ReadCountCollection factorNormalized = ReadCountCollectionUtils.parse(factorNormalizedOutput);
        final RealMatrix betaHats = readBetaHats(betaHatsOutput, input);

        Assert.assertEquals(factorNormalized.columnNames(), input.columnNames());
        Assert.assertTrue(!factorNormalized.targets().stream().anyMatch(t -> t.getInterval() == null));
        Assert.assertEquals(factorNormalized.targets().stream().map(Target::getInterval).collect(Collectors.toSet()),
               input.targets().stream().map(Target::getInterval).collect(Collectors.toSet()));
        Assert.assertEquals(factorNormalized.targets().stream().collect(Collectors.toSet()),
               input.targets().stream().collect(Collectors.toSet()));
        Assert.assertEquals(factorNormalized.columnNames(), input.columnNames());

        Assert.assertEquals(factorNormalized.targets().stream().collect(Collectors.toSet()),
               input.targets().stream().collect(Collectors.toSet()));
        assertFactorNormalizedValues(input, factorNormalized);

        final ReadCountCollection preTangentNormalized = ReadCountCollectionUtils.parse(preTangentNormalizationOutput);
        assertPreTangentNormalizedValues(factorNormalized, preTangentNormalized);
        assertBetaHats(preTangentNormalized, betaHats, TEST_PON);
        assertBetaHatsRobustToOutliers(preTangentNormalized, TEST_PON);

        // Test the tangent normalized output
        final ReadCountCollection tangentNormalized = ReadCountCollectionUtils.parse(tangentNormalizationOutput);
        assertTangentNormalized(tangentNormalized, preTangentNormalized, betaHats, TEST_PON);
        Assert.assertEquals(tangentNormalized.columnNames(), input.columnNames());
        Assert.assertEquals(tangentNormalized.columnNames(), preTangentNormalized.columnNames());

        // Make sure that we can read in the tangent normalized targets as a collection of TargetCoverage
        Assert.assertEquals(tangentNormalized.targets().size(), 657);
        Assert.assertEquals(tangentNormalized.targets().get(2).getName(), "target_179700_CRYBB1");
    }

    @Test(expectedExceptions = UserException.BadInput.class)
    public void testNameOnlyCountsInputRunMultipleSamples() throws IOException {
        final File factorNormalizedOutput = createTempFile("test",".txt");
        final File tangentNormalizationOutput = createTempFile("tangent-", ".txt");
        final File betaHatsOutput = createTempFile("tangent-", ".bhats");
        final File preTangentNormalizationOutput = createTempFile("pre-tn-",".txt");

        final String[] arguments = {
                "-" + NormalizeSomaticReadCounts.READ_COUNTS_FILE_SHORT_NAME, TARGET_NAME_ONLY_READ_COUNTS_INPUT.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.PON_FILE_SHORT_NAME, TEST_PON.getAbsolutePath(),
                "-" + NormalizeSomaticReadCounts.FACTOR_NORMALIZED_COUNTS_SHORT_NAME, factorNormalizedOutput.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.TANGENT_NORMALIZED_COUNTS_FILE_SHORT_NAME, tangentNormalizationOutput.getAbsolutePath(),
                "-" + NormalizeSomaticReadCounts.TANGENT_BETA_HATS_SHORT_NAME, betaHatsOutput.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.PRE_TANGENT_NORMALIZED_COUNTS_FILE_SHORT_NAME, preTangentNormalizationOutput.getAbsolutePath()
        };

        runCommandLine(arguments);
    }

    @Test
    public void testNameOnlyCountsInputRunOneSample() throws IOException {
        final File factorNormalizedOutput = createTempFile("test",".txt");
        final File tangentNormalizationOutput = createTempFile("tangent-", ".txt");
        final File betaHatsOutput = createTempFile("tangent-", ".bhats");
        final File preTangentNormalizationOutput = createTempFile("pre-tn-",".txt");

        final String[] arguments = {
                "-" + NormalizeSomaticReadCounts.READ_COUNTS_FILE_SHORT_NAME, TARGET_NAME_ONLY_READ_COUNTS_INPUT_ONE_SAMPLE.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.PON_FILE_SHORT_NAME, TEST_PON.getAbsolutePath(),
                "-" + NormalizeSomaticReadCounts.FACTOR_NORMALIZED_COUNTS_SHORT_NAME, factorNormalizedOutput.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.TANGENT_NORMALIZED_COUNTS_FILE_SHORT_NAME, tangentNormalizationOutput.getAbsolutePath(),
                "-" + NormalizeSomaticReadCounts.TANGENT_BETA_HATS_SHORT_NAME, betaHatsOutput.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.PRE_TANGENT_NORMALIZED_COUNTS_FILE_SHORT_NAME, preTangentNormalizationOutput.getAbsolutePath()
        };

        runCommandLine(arguments);
        final ReadCountCollection input = ReadCountCollectionUtils.parse(TARGET_NAME_ONLY_READ_COUNTS_INPUT_ONE_SAMPLE);
        final ReadCountCollection factorNormalized = ReadCountCollectionUtils.parse(factorNormalizedOutput);
        final ReadCountCollection tangentNormalized = ReadCountCollectionUtils.parse(tangentNormalizationOutput);
        final ReadCountCollection preTangentNormalized = ReadCountCollectionUtils.parse(preTangentNormalizationOutput);
        final RealMatrix betaHats = readBetaHats(betaHatsOutput, input);
        Assert.assertFalse(factorNormalized.targets().stream().anyMatch(t -> t.getInterval() != null));
        Assert.assertEquals(factorNormalized.columnNames(), input.columnNames());
        Assert.assertEquals(tangentNormalized.columnNames(), input.columnNames());
        Assert.assertEquals(preTangentNormalized.columnNames(), input.columnNames());
        Assert.assertEquals(factorNormalized.targets().stream().collect(Collectors.toSet()),
                input.targets().stream().collect(Collectors.toSet()));
        assertFactorNormalizedValues(input, factorNormalized);
        assertPreTangentNormalizedValues(factorNormalized, preTangentNormalized);
        assertBetaHats(preTangentNormalized, betaHats, TEST_PON);
        assertBetaHatsRobustToOutliers(preTangentNormalized, TEST_PON);
        assertTangentNormalized(tangentNormalized, preTangentNormalized, betaHats, TEST_PON);
    }

    private RealMatrix readBetaHats(final File betaHatsOutput, final ReadCountCollection input) throws IOException {
        final double[][] betaHats;

        try (final TableReader<double[]> reader = TableUtils.reader(betaHatsOutput,
                (columns, formatExceptionFactory) -> {
                    TableUtils.checkMandatoryColumns(columns, new TableColumnCollection(input.columnNames()), formatExceptionFactory);
                    if (!columns.matches(0, NormalizeSomaticReadCounts.PON_SAMPLE_BETA_HAT_COLUMN_NAME) ||
                            columns.columnCount() != input.columnNames().size() + 1) {
                        Assert.fail("Beta-hats has bad header");
                    }
                    return (dataLine) -> IntStream.range(0, input.columnNames().size())
                            .mapToDouble(i -> dataLine.getDouble(input.columnNames().get(i))).toArray();
                })) {
          betaHats = reader.stream().toArray(double[][]::new);
        }
        return new Array2DRowRealMatrix(betaHats,false);
    }

    @Test(expectedExceptions = UserException.BadInput.class)
    public void testCoordOnlyCountsRunMultipleSamples() throws IOException {
        final File factorNormalizedOutput = createTempFile("test",".txt");
        final File tangentNormalizationOutput = createTempFile("test",".txt");

        final String[] arguments = {
                "-" + NormalizeSomaticReadCounts.READ_COUNTS_FILE_SHORT_NAME, COORD_ONLY_READ_COUNTS_INPUT.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.PON_FILE_SHORT_NAME, TEST_PON.getAbsolutePath(),
                "-" + NormalizeSomaticReadCounts.FACTOR_NORMALIZED_COUNTS_SHORT_NAME, factorNormalizedOutput.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.TANGENT_NORMALIZED_COUNTS_FILE_SHORT_NAME, tangentNormalizationOutput.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.TARGET_FILE_SHORT_NAME, TEST_TARGETS.getAbsolutePath(),
        };

        runCommandLine(arguments);
    }

    @Test
    public void testCoordOnlyCountsInputFileRunOneSample() throws IOException {
        final File factorNormalizedOutput = createTempFile("test",".txt");
        final File tangentNormalizationOutput = createTempFile("test",".txt");

        final String[] arguments = {
                "-" + NormalizeSomaticReadCounts.READ_COUNTS_FILE_SHORT_NAME, COORD_ONLY_READ_COUNTS_INPUT_ONE_SAMPLE.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.PON_FILE_SHORT_NAME, TEST_PON.getAbsolutePath(),
                "-" + NormalizeSomaticReadCounts.FACTOR_NORMALIZED_COUNTS_SHORT_NAME, factorNormalizedOutput.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.TANGENT_NORMALIZED_COUNTS_FILE_SHORT_NAME, tangentNormalizationOutput.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.TARGET_FILE_SHORT_NAME, TEST_TARGETS.getAbsolutePath(),
        };

        runCommandLine(arguments);
        final TargetCollection<? extends BEDFeature> exons = TargetCollectionUtils.fromBEDFeatureFile(TEST_TARGETS, new BEDCodec());
        final ReadCountCollection input = ReadCountCollectionUtils.parse(COORD_ONLY_READ_COUNTS_INPUT_ONE_SAMPLE, exons, false);
        final ReadCountCollection factorNormalized = ReadCountCollectionUtils.parse(factorNormalizedOutput, exons, false);

        Assert.assertEquals(factorNormalized.columnNames(), input.columnNames());
        Assert.assertFalse(factorNormalized.targets().stream().anyMatch(t -> t.getInterval() == null));
       // Assert.assertTrue(factorNormalized.hasTargetNames());
       // Assert.assertEquals(factorNormalized.getIntervals(), input.getIntervals());
        assertFactorNormalizedValues(input, factorNormalized);
    }

    @Test(expectedExceptions = UserException.BadInput.class)
    public void testCoordOnlyCountsMissingTargetInputFileRunMultipleSamples() throws IOException {
        final File factorNormalizedOutput = createTempFile("test",".txt");
        final File tangentNormalizationOutput = createTempFile("test",".txt");

        final String[] arguments = {
                "-" + NormalizeSomaticReadCounts.READ_COUNTS_FILE_SHORT_NAME, COORD_ONLY_READ_COUNTS_INPUT.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.PON_FILE_SHORT_NAME, TEST_PON.getAbsolutePath(),
                "-" + NormalizeSomaticReadCounts.FACTOR_NORMALIZED_COUNTS_SHORT_NAME, factorNormalizedOutput.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.TANGENT_NORMALIZED_COUNTS_FILE_SHORT_NAME, tangentNormalizationOutput.getAbsolutePath(),
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
        try (final HDF5File ponReader = new HDF5File(TEST_PON)) {
            final PoN pon = new HDF5PoN(ponReader);
            final RealMatrix targetFactors = pon.getTargetFactors();
            final List<String> ponTargets = pon.getTargetNames();
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

    /**
     * This code reconstructs the expected pre-tangent normalization counts given the input,
     * and then compares against the actual pre-tangent output.
     * <p>
     * This method does not use any of the components that does the actual computation
     * in production, so that the calculation is independent.
     * </p>
     * <p>
     * This code also use an alternative way to calculate it to reduce overlap even further.
     * </p>
     * @param preTangentNormalized actual output.
     * @param factorNormalized input.
     */
    private void assertPreTangentNormalizedValues(final ReadCountCollection factorNormalized, final ReadCountCollection preTangentNormalized) {
        final double epsilon = TangentNormalizer.EPSILON;
        final RealMatrix outCounts = preTangentNormalized.counts();
        final RealMatrix inCounts = factorNormalized.counts();
        final double[] columnMeans = GATKProtectedMathUtils.columnMeans(inCounts);
        Assert.assertTrue(DoubleStream.of(columnMeans).allMatch(d -> d < 0.5));
        final double[][] expected = new double[inCounts.getRowDimension()][inCounts.getColumnDimension()];
        final double[] columnValues = new double[inCounts.getRowDimension()];
        for (int i = 0; i < columnMeans.length; i++) {
            for (int j = 0; j < inCounts.getRowDimension(); j++) {
                final double inValue = inCounts.getEntry(j,i);
                final double lowBoundedInValue = Math.max(epsilon, inValue / columnMeans[i]);
                final double outValue = Math.log(lowBoundedInValue) / Math.log(2);
                expected[j][i] = outValue;
                columnValues[j] = outValue;
            }
            Arrays.sort(columnValues);
            final int midIndex = columnValues.length >> 1;
            final double median = columnValues.length % 2 == 1 ?
                    columnValues[midIndex] : (columnValues[midIndex] + columnValues[1 + midIndex]) * .5;
            for (int j = 0; j < inCounts.getRowDimension(); j++) {
                expected[j][i] -= median;
                Assert.assertEquals(outCounts.getEntry(j,i),expected[j][i],0.000001," Row " + j + " col " + i);
            }
        }
    }

    /**
     * Asserts that the calculation of beta hats is not significantly affected by zero-coverage outlier counts
     * We perform this check by randomly setting some coverages to zero in copy ratio space (-infinity in log space).
     * betaHats imputes 0 in log space (1 in copy ratio space) whenever coverage is below a certain low threshold
     * and should thus be robust to this type of noise.
     */
    private void assertBetaHatsRobustToOutliers(final ReadCountCollection preTangentNormalized, final File ponFile) {
        try (final HDF5File ponReader = new HDF5File(ponFile)) {
            final PoN pon = new HDF5PoN(ponReader);
            final List<String> ponTargets = pon.getPanelTargetNames();

            final RealMatrix input = reorderTargetsToPoNOrder(preTangentNormalized, ponTargets);

            // randomly set some entries to zero in copy-ratio space (-infinity in log space)
            final Random random = new Random(13);
            final double noiseProportion = 0.01;
            final RealMatrix noisyInput = input.copy();
            noisyInput.walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor() {
                @Override
                public double visit(final int row, final int column, final double value) {
                    return random.nextDouble() < noiseProportion ? Double.NEGATIVE_INFINITY : value;
                }
            });

            final RealMatrix betaHats = pon.betaHats(input, true, TangentNormalizer.EPSILON);
            final RealMatrix noisyBetaHats = pon.betaHats(noisyInput, true, TangentNormalizer.EPSILON);
            final RealMatrix difference = betaHats.subtract(noisyBetaHats);

            difference.walkInOptimizedOrder(new DefaultRealMatrixPreservingVisitor() {
                @Override
                public void visit(final int row, int column, double value) {
                    Assert.assertEquals(value, 0, 0.01);
                }
            });
        }
    }

    /**
     * Asserts that a collection of beta-hats corresponds to the expected value given
     * the input pre-tangent normalization matrix and the PoN file.
     */
    private void assertBetaHats(final ReadCountCollection preTangentNormalized,
                                final RealMatrix actual, final File ponFile) {
        Assert.assertEquals(actual.getColumnDimension(), preTangentNormalized.columnNames().size());
        final double epsilon = TangentNormalizer.EPSILON;

        try (final HDF5File ponReader = new HDF5File(ponFile)) {
            final PoN pon = new HDF5PoN(ponReader);
            final List<String> ponTargets = pon.getPanelTargetNames();
            final RealMatrix inCounts = reorderTargetsToPoNOrder(preTangentNormalized, ponTargets);

            // obtain subset of relevant targets to calculate the beta-hats;
            final int[][] betaHatTargets = new int[inCounts.getColumnDimension()][];
            for (int i = 0; i < inCounts.getColumnDimension(); i++) {
                final List<Integer> relevantTargets = new ArrayList<>();
                for (int j = 0; j < inCounts.getRowDimension(); j++) {
                    if (inCounts.getEntry(j, i) > 1  +  (Math.log(epsilon) / Math.log(2))) {
                        relevantTargets.add(j);
                    }
                }
                betaHatTargets[i] = relevantTargets.stream().mapToInt(Integer::intValue).toArray();
            }
            // calculate beta-hats per column and check with actual values.
            final RealMatrix normalsInv = pon.getReducedPanelPInverseCounts();
            Assert.assertEquals(actual.getRowDimension(), normalsInv.getRowDimension());
            final RealMatrix normalsInvT = normalsInv.transpose();
            for (int i = 0; i < inCounts.getColumnDimension(); i++) {
                final RealMatrix inValues = inCounts.getColumnMatrix(i).transpose().getSubMatrix(new int[] { 0 }, betaHatTargets[i]);
                final RealMatrix normalValues = normalsInvT.getSubMatrix(betaHatTargets[i], IntStream.range(0, normalsInvT.getColumnDimension()).toArray());
                final RealMatrix betaHats = inValues.multiply(normalValues);
                for (int j = 0; j < actual.getRowDimension(); j++) {
                    Assert.assertEquals(actual.getEntry(j, i), betaHats.getEntry(0, j),0.000001,"Col " + i + " row " + j);
                }
            }
        }
    }

    private void assertTangentNormalized(final ReadCountCollection actualReadCounts, final ReadCountCollection preTangentNormalized, final RealMatrix betaHats, final File ponFile) {

        try (final HDF5File ponReader = new HDF5File(ponFile)) {
            final PoN pon = new HDF5PoN(ponReader);
            final RealMatrix inCounts = reorderTargetsToPoNOrder(preTangentNormalized, pon.getPanelTargetNames());
            final RealMatrix actual = reorderTargetsToPoNOrder(actualReadCounts,pon.getPanelTargetNames());
            final RealMatrix ponMat = pon.getReducedPanelCounts();
            final RealMatrix projection = ponMat.multiply(betaHats);
            final RealMatrix expected = inCounts.subtract(projection);
            Assert.assertEquals(actual.getRowDimension(),expected.getRowDimension());
            Assert.assertEquals(actual.getColumnDimension(),expected.getColumnDimension());
            for (int i = 0; i < actual.getRowDimension(); i++) {
                Assert.assertEquals(actual.getRow(i),expected.getRow(i));
            }
        }
    }

    private RealMatrix reorderTargetsToPoNOrder(final ReadCountCollection preTangentNormalized, final List<String> ponTargets) {
        final RealMatrix preTangentNormalizedCounts = preTangentNormalized.counts();
        final Map<String,Integer> ponTargetIndex = IntStream.range(0, ponTargets.size())
                .boxed().collect(Collectors.toMap(ponTargets::get, Function.identity()));

        // first we need to sort the input counts so that they match the
        // target order in the PoN.
        final double[][] ponPreparedInput = new double[ponTargets.size()][];
        for (int i = 0; i < preTangentNormalizedCounts.getRowDimension(); i++) {
            final Target target = preTangentNormalized.targets().get(i);
            if (!ponTargetIndex.containsKey(target.getName()))
                continue;
            final int idx = ponTargetIndex.get(target.getName());
            ponPreparedInput[idx] = preTangentNormalizedCounts.getRow(i);
        }

        // The actual input to create the beta-hats, sorted by the PoN targets:
        return new Array2DRowRealMatrix(ponPreparedInput,false);
    }
}
