package org.broadinstitute.hellbender.tools.exome;

import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.ArgumentCollection;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SparkToggleCommandLineProgram;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.hdf5.*;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.io.File;
import java.util.*;

/**
 * Tool to create a panel of normals (PoN) given a collection of read-counts
 * for control samples.
 *
 * <p>
 * The input read-counts consists of a single file with counts for several
 * samples. This might be constructed from several single sample read counts using
 * the {@link CombineReadCounts} tool.
 * </p>
 * <p>
 * Accordingly the input format is exactly the same as the output format for
 * that tool, which is described in its documentation {@link CombineReadCounts here}.
 * </p>
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "Creates a Panel of Normals (PoN) given the proportional read counts for the samples that are part of the panel.  Supports Apache Spark for some operations.",
        oneLineSummary = "Creates a Panel of Normals.",
        programGroup = CopyNumberProgramGroup.class
)
public class CreatePanelOfNormals extends SparkToggleCommandLineProgram {

    static final long serialVersionUID = 42123132L;

    public static final double DEFAULT_TARGET_FACTOR_THRESHOLD_PERCENTILE = 25.0;

    public static final double DEFAULT_COLUMN_OUTLIER_DROP_THRESHOLD_PERCENTILE = 2.5;

    // 0.1% by default as copied from python code... perhaps a bug there as most of these constants in that code
    // are then multiplied by 100 or subtracted from a 100, check out PanelCleaner.py in ReCapSeg repo.
    public static final double DEFAULT_OUTLIER_TRUNCATE_PERCENTILE_THRESHOLD = 0.1;

    public static final double DEFAULT_MAXIMUM_PERCENT_ZEROS_IN_COLUMN = 2.0;

    public static final double DEFAULT_MAXIMUM_PERCENT_ZEROS_IN_TARGET = 5.0;

    public static final String INFER_NUMBER_OF_EIGEN_SAMPLES = "auto";

    public static final String DEFAULT_NUMBER_OF_EIGEN_SAMPLES = INFER_NUMBER_OF_EIGEN_SAMPLES;


    public static final String NUMBER_OF_EIGEN_SAMPLES_DOCUMENTATION =
            "Number of eigen samples to use for the reduced PoN. " +
            "By default it will infer the appropriate number of eigen samples (value " + INFER_NUMBER_OF_EIGEN_SAMPLES + ")";

    public static final String COLUMN_EXTREME_THRESHOLD_PERCENTILE_DOCUMENTATION =
            "Percentile for the two-tailed extreme median column coverage filter. " +
            "Columns that have a median count in the that bottom or top percentile are excluded from the PoN. " +
            "Any value in the range (0, 50) (default is " + DEFAULT_COLUMN_OUTLIER_DROP_THRESHOLD_PERCENTILE + ")";

    public static final String TARGET_FACTOR_THRESHOLD_PERCENTILE_DOCUMENTATION =
            "Percentile to determine the minimum target factor for any target to be considered part of the PoN. " +
            "Any value in the range (0, 100) (default is " + DEFAULT_TARGET_FACTOR_THRESHOLD_PERCENTILE + ")";

    public static final String MAXIMUM_PERCENT_ZEROS_IN_COLUMN_DOCUMENTATION =
            "Maximum percentage of 0 counts in a column to drop it from the panel. " +
            "Any value in the range (0, 100) (default is " + DEFAULT_MAXIMUM_PERCENT_ZEROS_IN_COLUMN + ")";

    public static final String MAXIMUM_PERCENT_ZEROS_IN_TARGET_DOCUMENTATION =
            "Maximum percentage of 0 counts for a target to drop it from the panel. " +
                    "Any value in the range (0, 100) (default is " + DEFAULT_MAXIMUM_PERCENT_ZEROS_IN_TARGET + ")";

    public static final String COUNT_TRUNCATE_PERCENTILE_DOCUMENTATION =
            "Percentiles to obtain the maximum and minimum value for any count in the panel. " +
            "Any value outside the resulting range would be set will be truncated to that minimum or maximum. " +
            "Valid values are within the range (0, 50) (default is " + DEFAULT_OUTLIER_TRUNCATE_PERCENTILE_THRESHOLD + ")";


    public static final String TARGET_FACTOR_THRESHOLD_PERCENTILE_SHORT_NAME = "minTFPcTh";
    public static final String TARGET_FACTOR_THRESHOLD_PERCENTILE_FULL_NAME = "minimumTargetFactorPercentileThreshold";
    public static final String MAXIMUM_PERCENT_ZEROS_IN_COLUMN_SHORT_NAME = "maxCol0sPc";
    public static final String MAXIMUM_PERCENT_ZEROS_IN_COLUMN_FULL_NAME = "maximumColumnZerosPercentage";
    public static final String MAXIMUM_PERCENT_ZEROS_IN_TARGET_SHORT_NAME = "maxTrg0sPc";
    public static final String MAXIMUM_PERCENT_ZEROS_IN_TARGET_FULL_NAME = "maximumTargetZerosPercentage";
    public static final String COLUMN_EXTREME_THRESHOLD_PERCENTILE_SHORT_NAME = "extremeColMedPrTh";
    public static final String COLUMN_EXTREME_THRESHOLD_PERCENTILE_FULL_NAME = "extremeColumnMedianCountPercentileThreshold";
    public static final String COUNT_TRUNCATE_PERCENTILE_SHORT_NAME = "truncPcThr";
    public static final String COUNT_TRUNCATE_PERCENTILE_FULL_NAME = "truncatePercentileThreshold";
    public static final String NUMBER_OF_EIGEN_SAMPLES_SHORT_NAME = "numEigen";
    public static final String NUMBER_OF_EIGEN_SAMPLES_FULL_NAME = "numberOfEigenSamples";
    public static final String DRY_RUN_SHORT_NAME = "dryRun";
    public static final String DRY_RUN_FULL_NAME = DRY_RUN_SHORT_NAME;
    public static final String NO_QC_SHORT_NAME = "noQC";
    public static final String NO_QC_FULL_NAME = NO_QC_SHORT_NAME;
    public static final String TARGET_WEIGHTS_SHORT_NAME = "tw";
    public static final String TARGET_WEIGHTS_FULL_NAME = "outputTargetWeights";
    public static final String BLACKLIST_QC_SHORT_NAME = "bo";
    public static final String BLACKLIST_QC_FULL_NAME = "outputFailedSamples";
    public static final int NUM_QC_EIGENSAMPLES = 20;
    public static final String BLACKLIST_FILE_APPEND = ".removed_samples.txt";
    public static final String TARGET_WEIGHTS_FILE_APPEND = ".target_weights.txt";


    @Argument(
            doc = TARGET_FACTOR_THRESHOLD_PERCENTILE_DOCUMENTATION,
            shortName = TARGET_FACTOR_THRESHOLD_PERCENTILE_SHORT_NAME,
            fullName  = TARGET_FACTOR_THRESHOLD_PERCENTILE_FULL_NAME,
            optional  = true
    )
    protected double targetFactorThreshold = DEFAULT_TARGET_FACTOR_THRESHOLD_PERCENTILE;

    @Argument(
            doc = MAXIMUM_PERCENT_ZEROS_IN_COLUMN_DOCUMENTATION,
            shortName = MAXIMUM_PERCENT_ZEROS_IN_COLUMN_SHORT_NAME,
            fullName = MAXIMUM_PERCENT_ZEROS_IN_COLUMN_FULL_NAME,
            optional = true
    )
    protected double maximumPercentZerosInColumn = DEFAULT_MAXIMUM_PERCENT_ZEROS_IN_COLUMN;

    @Argument(
            doc = MAXIMUM_PERCENT_ZEROS_IN_TARGET_DOCUMENTATION,
            shortName = MAXIMUM_PERCENT_ZEROS_IN_TARGET_SHORT_NAME,
            fullName = MAXIMUM_PERCENT_ZEROS_IN_TARGET_FULL_NAME,
            optional = true
    )
    protected double maximumPercentZerosInTarget = DEFAULT_MAXIMUM_PERCENT_ZEROS_IN_TARGET;

    @Argument(
            doc = COLUMN_EXTREME_THRESHOLD_PERCENTILE_DOCUMENTATION,
            shortName = COLUMN_EXTREME_THRESHOLD_PERCENTILE_SHORT_NAME,
            fullName = COLUMN_EXTREME_THRESHOLD_PERCENTILE_FULL_NAME,
            optional = false
    )
    protected double columnExtremeThresholdPercentile = DEFAULT_COLUMN_OUTLIER_DROP_THRESHOLD_PERCENTILE;

    @Argument(
            doc = COUNT_TRUNCATE_PERCENTILE_DOCUMENTATION,
            shortName = COUNT_TRUNCATE_PERCENTILE_SHORT_NAME,
            fullName = COUNT_TRUNCATE_PERCENTILE_FULL_NAME,
            optional = false
    )
    protected double outlierTruncatePercentileThresh = DEFAULT_OUTLIER_TRUNCATE_PERCENTILE_THRESHOLD;

    @Argument(
            doc = NUMBER_OF_EIGEN_SAMPLES_DOCUMENTATION,
            shortName = NUMBER_OF_EIGEN_SAMPLES_SHORT_NAME,
            fullName = NUMBER_OF_EIGEN_SAMPLES_FULL_NAME,
            optional = true
    )
    protected String numberOfEigenSamples = DEFAULT_NUMBER_OF_EIGEN_SAMPLES;

    @Argument(
            doc = "Skip the QC step.  PoN creation will be substantially faster, but greater risk of bad samples being introduced into the PoN.",
            shortName = NO_QC_SHORT_NAME,
            fullName  = NO_QC_FULL_NAME,
            optional  = true
    )
    protected boolean isNoQc = false;

    @Argument(
            doc = "Input proportional read counts for samples in the panel of normals.",
            shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME,
            fullName  = StandardArgumentDefinitions.INPUT_LONG_NAME,
            optional  = false
    )
    protected File inputFile = null;

    @ArgumentCollection
    protected TargetArgumentCollection targetArguments = new TargetArgumentCollection(() -> inputFile);

    @Argument(
            doc = "Output HDF5 file name.",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            optional = false
    )
    protected File outFile = null;

    @Argument(
            doc = "Output file for sample names that failed the QC check (one sample name per line; sample names drawn from input file).  This option is ignored if QC step is being skipped.  Final PoN will not include these samples, if QC is enabled.",
            shortName = BLACKLIST_QC_SHORT_NAME,
            fullName = BLACKLIST_QC_FULL_NAME,
            optional = true
    )
    protected File blacklistOutFile = null;

    @Argument(
            doc = "Output file for the weights of each target.  This is the 1/var(post-projected targets for each normal)." +
                    "  This file can be used for the PerformSegmentationTool if desired.  By default, this simple text file is " +
                    "written next to the pon file with the extension " + TARGET_WEIGHTS_FILE_APPEND + ".",
            shortName = TARGET_WEIGHTS_SHORT_NAME,
            fullName = TARGET_WEIGHTS_FULL_NAME,
            optional = true
    )
    protected File targetWeightsOutFile = null;

    // This option is useful to test performance when the HDF5 lib is not
    // present for whatever reason.
    @Argument(
            doc = "Dry-run, skip the creation of the HDF5 output file.  Will also skip QC, if specified.",
            shortName = DRY_RUN_SHORT_NAME,
            fullName  = DRY_RUN_FULL_NAME,
            optional  = true
    )
    protected boolean dryRun = false;

    @Override
    protected void runPipeline(JavaSparkContext ctx) {

        if (blacklistOutFile == null) {
            blacklistOutFile = new File(outFile + BLACKLIST_FILE_APPEND);
        }
        if (targetWeightsOutFile == null) {
            targetWeightsOutFile = new File(outFile + TARGET_WEIGHTS_FILE_APPEND);
        }

        // Check parameters and load values to meet the backend PoN creation interface
        final TargetCollection<Target> targets = targetArguments.readTargetCollection(true);
        final double targetFactorPercentileThreshold = checkTargetFactorsPercentileThreshold();
        final double extremeColumnMedianCountPercentileThreshold = checkExtremeColumnMedianCountsPercentileThreshold();
        final double countTruncatePercentile = checkCountTruncatePercentile();
        final double maximumPercentageZeroTargets = checkTargetMaximumZeroPercentage();
        final double maximumPercentageZeroColumns = checkColumnMaximumZeroPercentage();
        final OptionalInt numberOfEigenSamples = calculatePreferredNumberOfEigenSamples();

        // Create the PoN, including QC, if specified.
        if (!isNoQc && !dryRun) {
            logger.info("QC:  Beginning creation of QC PoN...");
            final File outputQCFile = IOUtils.createTempFile("qc-pon-",".hd5");
            HDF5PoNCreator.createPoN(ctx, inputFile, OptionalInt.of(NUM_QC_EIGENSAMPLES), new ArrayList<>(),
                    outputQCFile, dryRun, targets, targetFactorPercentileThreshold, extremeColumnMedianCountPercentileThreshold,
                    countTruncatePercentile, maximumPercentageZeroTargets, maximumPercentageZeroColumns);
            logger.info("QC:  QC PoN created...");

            logger.info("QC:  Collecting suspicious samples...");
            final List<String> failingSampleNames = PoNIssueDetector.retrieveSamplesWithArmLevelEvents(outputQCFile, ctx);

            ParamUtils.writeStringListToFile(failingSampleNames, blacklistOutFile);

            // If no suspicious samples were found, just redo the PoN reduction to save time.
            if (failingSampleNames.size() != 0) {
                logger.info("QC:  Suspicious sample list created...");

                logger.info("Creating final PoN with " + failingSampleNames.size() + " suspicious samples removed...");
                HDF5PoNCreator.createPoN(ctx, inputFile, numberOfEigenSamples, failingSampleNames, outFile, dryRun,
                        targets, targetFactorPercentileThreshold, extremeColumnMedianCountPercentileThreshold,
                        countTruncatePercentile, maximumPercentageZeroTargets, maximumPercentageZeroColumns);
            } else {
                logger.info("QC:  No suspicious samples found ...");
                logger.info("Creating final PoN only redo'ing the reduction step ...");
                HDF5PoNCreator.redoReduction(ctx, numberOfEigenSamples, outputQCFile, outFile);
            }

        } else {
            logger.info("Creating PoN directly (skipping QC)...");
            HDF5PoNCreator.createPoN(ctx, inputFile, numberOfEigenSamples, new ArrayList<>(), outFile, dryRun, targets,
                    targetFactorPercentileThreshold, extremeColumnMedianCountPercentileThreshold,
                    countTruncatePercentile, maximumPercentageZeroTargets, maximumPercentageZeroColumns);
        }

        if (!dryRun) {
            logger.info("Writing target weights file to " + targetWeightsOutFile + "...");
            writeTargetWeightsFile(outFile, targetWeightsOutFile);
        }
        logger.info("Done...");
    }


    /**
     * Composes the preferred number of eigen values optional given the user input.
     *
     * @return an empty optional if the user elected to use the automatic/inferred value, otherwise
     * a strictly positive integer.
     */
    private OptionalInt calculatePreferredNumberOfEigenSamples() {
        if (numberOfEigenSamples.equalsIgnoreCase(INFER_NUMBER_OF_EIGEN_SAMPLES)) {
            return OptionalInt.empty();
        } else {
            try {
                final int result = Integer.parseInt(numberOfEigenSamples);
                if (result <= 0) {
                    throw new UserException.BadArgumentValue(NUMBER_OF_EIGEN_SAMPLES_FULL_NAME, "0 or negative values are not allowed: " + numberOfEigenSamples);
                } else {
                    return OptionalInt.of(result);
                }
            } catch (final NumberFormatException ex) {
                throw new UserException.BadArgumentValue(NUMBER_OF_EIGEN_SAMPLES_FULL_NAME,
                        "it must be either '" + INFER_NUMBER_OF_EIGEN_SAMPLES + "' or an integer value");
            }
        }
    }

    private double checkExtremeColumnMedianCountsPercentileThreshold() {
        if (columnExtremeThresholdPercentile < 0 || columnExtremeThresholdPercentile > 50 || Double.isNaN(columnExtremeThresholdPercentile)) {
            throw new UserException.BadArgumentValue(COLUMN_EXTREME_THRESHOLD_PERCENTILE_FULL_NAME, "the value must be in the range [0, 50]");
        }
        return columnExtremeThresholdPercentile;
    }

    private double checkTargetFactorsPercentileThreshold() {
        if (targetFactorThreshold < 0 || targetFactorThreshold > 100 || Double.isNaN(targetFactorThreshold)) {
            throw new UserException.BadArgumentValue(TARGET_FACTOR_THRESHOLD_PERCENTILE_FULL_NAME, "the value must be in the range [0, 100]");
        }
        return targetFactorThreshold;
    }

    /**
     * Checks that the user input {@link #outlierTruncatePercentileThresh} is correct.
     * @return the user input {@link #outlierTruncatePercentileThresh}.
     */
    private double checkCountTruncatePercentile() {
        if (outlierTruncatePercentileThresh < 0 || outlierTruncatePercentileThresh > 50.0 || Double.isNaN(outlierTruncatePercentileThresh)) {
            throw new UserException.BadArgumentValue(COUNT_TRUNCATE_PERCENTILE_FULL_NAME, "the value must be in the [0, 50.0) range");
        }
        return outlierTruncatePercentileThresh;
    }

    private double checkTargetMaximumZeroPercentage() {
        if (maximumPercentZerosInTarget < 0 || maximumPercentZerosInTarget > 100 || Double.isNaN(maximumPercentZerosInTarget)) {
            throw new UserException.BadArgumentValue(MAXIMUM_PERCENT_ZEROS_IN_TARGET_FULL_NAME, "the value must be in the [0, 100.0) range");
        }
        return maximumPercentZerosInTarget;
    }

    private double checkColumnMaximumZeroPercentage() {
        if (maximumPercentZerosInColumn < 0 || maximumPercentZerosInColumn > 100 || Double.isNaN(maximumPercentZerosInColumn)) {
            throw new UserException.BadArgumentValue(MAXIMUM_PERCENT_ZEROS_IN_COLUMN_FULL_NAME, "the value must be in the [0, 100.0) range");
        }
        return maximumPercentZerosInColumn;
    }

    /**
     * Write a target weight file that can be read in by the R CBS.
     *
     * @param pon never {@code null}
     * @param outputFile never {@code null}.  Must be writable.
     */
    public static void writeTargetWeightsFile(final PoN pon, final File outputFile) {
        final double[] targetWeights = HDF5PoNCreator.calculateTargetWeights(pon);
        ParamUtils.writeValuesToFile(targetWeights, outputFile);
    }

    /**
     * Write a target weight file that can be read in by the R CBS.
     *
     * @param ponFile never {@code null}.  Must be readable.
     * @param outputFile never {@code null}.  Must be writable.
     */
    public static void writeTargetWeightsFile(final File ponFile, final File outputFile) {
        Utils.regularReadableUserFile(ponFile);
        try (final HDF5File file = new HDF5File(ponFile, HDF5File.OpenMode.READ_ONLY)) {
            final HDF5PoN pon = new HDF5PoN(file);
            writeTargetWeightsFile(pon, outputFile);
        }
    }
}
