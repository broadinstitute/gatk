package org.broadinstitute.hellbender.tools.walkers.vqsr.scalable;

import com.google.common.collect.Streams;
import com.google.common.primitives.Doubles;
import org.apache.commons.math3.stat.descriptive.moment.Variance;
import org.apache.commons.math3.stat.descriptive.rank.Percentile;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.vqsr.VariantRecalibrator;
import org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.data.LabeledVariantAnnotationsData;
import org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.data.VariantType;
import org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.modeling.BGMMVariantAnnotationsModel;
import org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.modeling.BGMMVariantAnnotationsScorer;
import org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.modeling.PythonSklearnVariantAnnotationsModel;
import org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.modeling.PythonSklearnVariantAnnotationsScorer;
import org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.modeling.VariantAnnotationsModel;
import org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.modeling.VariantAnnotationsModelBackend;
import org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.modeling.VariantAnnotationsScorer;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.io.Resource;
import org.broadinstitute.hellbender.utils.python.PythonScriptExecutor;
import picard.cmdline.programgroups.VariantFilteringProgramGroup;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Trains a model for scoring variant calls based on site-level annotations.
 *
 * <p>
 *     This tool is intended to be used as the second step in a variant-filtering workflow that supersedes the
 *     {@link VariantRecalibrator} workflow. Given training (and optionally, calibration) sets of site-level annotations
 *     produced by {@link ExtractVariantAnnotations}, this tool can be used to train a model for scoring variant
 *     calls. The outputs of the tool are TODO
 * </p>
 *
 * <p>
 *     The model trained by this tool can in turn be provided along with a VCF file to the {@link ScoreVariantAnnotations}
 *     tool, which assigns a score to each call (with a lower score indicating that a call is more likely to be an artifact
 *     and should perhaps be filtered). Each score can also be converted to a corresponding sensitivity to a
 *     calibration set, if the latter is available.
 * </p>
 *
 * <p>
 *     TODO model definition
 * </p>
 *
 * <p>
 *     TODO calibration-sensitivity conversion, considerations, and comparison to tranche files
 * </p>
 *
 * <p>
 *     TODO positive vs. positive-negative
 * </p>
 *  *
 * <p>
 *     TODO IsolationForest section with description of method and hyperparameters
 * </p>
 *
 * <p>
 *     Note that HDF5 files may be viewed using <a href="https://support.hdfgroup.org/products/java/hdfview/">hdfview</a>
 *     or loaded in Python using <a href="http://www.pytables.org/">PyTables</a> or <a href="http://www.h5py.org/">h5py</a>.
 * </p>
 *
 * <h3>Inputs</h3>
 *
 * <ul>
 *     <li>
 *         Labeled-annotations HDF5 file (.annot.hdf5). Annotation data and metadata for labeled sites are stored in the
 *         HDF5 directory structure given in the documentation for the {@link ExtractVariantAnnotations} tool. In typical
 *         usage, both the {@value LabeledVariantAnnotationsData#TRAINING_LABEL} and
 *         {@value LabeledVariantAnnotationsData#CALIBRATION_LABEL} labels would be available for non-empty sets of
 *         sites of the requested variant type.
 *     </li>
 *     <li>
 *         (Optional) Unlabeled-annotations HDF5 file (.unlabeled.annot.hdf5). Annotation data and metadata for
 *         unlabeled sites are stored in the HDF5 directory structure given in the documentation for the
 *         {@link ExtractVariantAnnotations} tool. If provided, a positive-negative modeling approach (similar to
 *         that used in {@link VariantRecalibrator} will be used.
 *     </li>
 *     <li>
 *         Variant types (i.e., SNP and/or INDEL) for which to train models. Logic for determining variant type was retained from
 *         {@link VariantRecalibrator}; see {@link VariantType}. A separate model will be trained for each variant type
 *         and separate sets of outputs with corresponding tags in the filenames (i.e., "snp" or "indel") will be produced.
 *         TODO can run tool twice
 *     </li>
 *     <li>
 *         (Optional) Model backend. The default Python IsolationForest implementation requires either the GATK Python environment
 *         or that certain Python packages (argparse, h5py, numpy, sklearn, and dill) are otherwise available.
 *         A custom backend can also be specified in conjunction with the {@value PYTHON_SCRIPT_LONG_NAME} argument.
 *     </li>
 *     <li>
 *         (Optional) Model hyperparameters JSON file. TODO
 *     </li>
 *     <li>
 *         (Optional) Calibration-set sensitivity threshold. TODO if separate SNP/INDEL thresholds, run tool twice
 *     </li>
 *     <li>
 *         Output prefix.
 *         This is used as the basename for output files.
 *     </li>
 * </ul>
 *
 * <h3>Outputs</h3>
 *
 * <ul>
 *     <li>
 *         TODO
 *     </li>
 *     <li>
 *         (Optional) TODO
 *     </li>
 * </ul>
 *
 * <h3>Usage examples</h3>
 *
 * <p>
 *     TODO, positive-only, producing the outputs 1)
 *
 * <pre>
 *     gatk TrainVariantAnnotationsModel \
 *          TODO
 * </pre>
 * </p>
 *
 * <p>
 *     TODO, positive-negative, producing the outputs 1)
 *
 * <pre>
 *     gatk TrainVariantAnnotationsModel \
 *          TODO
 * </pre>
 * </p>
 *
 * <h3>Custom modeling/scoring backends (ADVANCED)</h3>
 *
 * <p>
 *     The primary modeling functionality performed by this tool is accomplished by a "modeling backend"
 *     whose fundamental contract is to take an input HDF5 file containing an annotation matrix for sites of a
 *     single variant type (i.e., SNP or INDEL) and to output a serialized scorer for that variant type.
 *     Rather than using one of the available, implemented backends, advanced users may provide their own backend
 *     via the {@value PYTHON_SCRIPT_LONG_NAME} argument. See documentation in the modeling and scoring interfaces
 *     ({@link VariantAnnotationsModel} and {@link VariantAnnotationsScorer}, respectively), as well as the default
 *     Python IsolationForest implementation at {@link PythonSklearnVariantAnnotationsModel} and
 *     org/broadinstitute/hellbender/tools/walkers/vqsr/scalable/isolation-forest.py.
 * </p>
 *
 * <p>
 *     Extremely advanced users could potentially substitute their own implementation for the entire
 *     {@link TrainVariantAnnotationsModel} tool, while still making use of the up/downstream
 *     {@link ExtractVariantAnnotations} and {@link ScoreVariantAnnotations} tools. To do so, one would additionally
 *     have to implement functionality for subsetting training/calibration sets by variant type,
 *     calling modeling backends as appropriate, and scoring calibration sets.
 * </p>
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "Trains a model for scoring variant calls based on site-level annotations.",
        oneLineSummary = "Trains a model for scoring variant calls based on site-level annotations",
        programGroup = VariantFilteringProgramGroup.class
)
@DocumentedFeature
@BetaFeature
public final class TrainVariantAnnotationsModel extends CommandLineProgram {

    public static final String MODE_LONG_NAME = "mode";
    public static final String ANNOTATIONS_HDF5_LONG_NAME = "annotations-hdf5";
    public static final String UNLABELED_ANNOTATIONS_HDF5_LONG_NAME = "unlabeled-annotations-hdf5";
    public static final String MODEL_BACKEND_LONG_NAME = "model-backend";
    public static final String PYTHON_SCRIPT_LONG_NAME = "python-script";
    public static final String HYPERPARAMETERS_JSON_LONG_NAME = "hyperparameters-json";
    public static final String CALIBRATION_SENSITIVITY_THRESHOLD_LONG_NAME = "calibration-sensitivity-threshold";

    public static final String ISOLATION_FOREST_PYTHON_SCRIPT = "isolation-forest.py";
    public static final String ISOLATION_FOREST_HYPERPARAMETERS_JSON = "isolation-forest-hyperparameters.json";

    enum AvailableLabelsMode {
        POSITIVE_ONLY, POSITIVE_UNLABELED
    }

    public static final String TRAINING_SCORES_HDF5_SUFFIX = ".trainingScores.hdf5";
    public static final String CALIBRATION_SCORES_HDF5_SUFFIX = ".calibrationScores.hdf5";
    public static final String UNLABELED_SCORES_HDF5_SUFFIX = ".unlabeledScores.hdf5";
    public static final String NEGATIVE_TAG = ".negative";

    @Argument(
            fullName = ANNOTATIONS_HDF5_LONG_NAME,
            doc = "HDF5 file containing annotations extracted with ExtractVariantAnnotations.")
    private File inputAnnotationsFile;

    @Argument(
            fullName = UNLABELED_ANNOTATIONS_HDF5_LONG_NAME,
            doc = "HDF5 file containing annotations extracted with ExtractVariantAnnotations. " +
                    "If specified with " + CALIBRATION_SENSITIVITY_THRESHOLD_LONG_NAME + ", " +
                    "a positive-unlabeled modeling approach will be used; otherwise, a positive-only modeling " +
                    "approach will be used.",
            optional = true)
    private File inputUnlabeledAnnotationsFile;

    @Argument(
            fullName = MODEL_BACKEND_LONG_NAME,
            doc = "Backend to use for training models. " +
                    "JAVA_BGMM will use a pure Java implementation (ported from Python scikit-learn) of the Bayesian Gaussian Mixture Model. " +
                    "PYTHON_IFOREST will use the Python scikit-learn implementation of the IsolationForest method and " +
                    "will require that the corresponding Python dependencies are present in the environment. " +
                    "PYTHON_SCRIPT will use the script specified by the " + PYTHON_SCRIPT_LONG_NAME + " argument. " +
                    "See the tool documentation for more details.")
    private VariantAnnotationsModelBackend modelBackend = VariantAnnotationsModelBackend.PYTHON_IFOREST;

    @Argument(
            fullName = PYTHON_SCRIPT_LONG_NAME,
            doc = "Python script used for specifying a custom scoring backend. If provided, " + MODEL_BACKEND_LONG_NAME + " must also be set to PYTHON_SCRIPT.",
            optional = true)
    private File pythonScriptFile;

    @Argument(
            fullName = HYPERPARAMETERS_JSON_LONG_NAME,
            doc = "JSON file containing hyperparameters. Optional if the PYTHON_IFOREST backend is used " +
                    "(if not specified, a default set of hyperparameters will be used); otherwise required.",
            optional = true)
    private File hyperparametersJSONFile;

    @Argument(
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "Output prefix.")
    private String outputPrefix;

    @Argument(
            fullName = CALIBRATION_SENSITIVITY_THRESHOLD_LONG_NAME,
            doc = "Calibration-sensitivity threshold that determines which sites will be used for training the negative model " +
                    "in the positive-unlabeled modeling approach. " +
                    "Increasing this will decrease the corresponding positive-model score threshold; sites with scores below this score " +
                    "threshold will be used for training the negative model. Thus, this parameter should typically be chosen to " +
                    "be close to 1, so that sites that score highly according to the positive model will not be used to train the negative model. " +
                    "The " + UNLABELED_ANNOTATIONS_HDF5_LONG_NAME + " argument must be specified in conjunction with this argument. " +
                    "If separate thresholds for SNP and INDEL models are desired, run the tool separately for each mode with its respective threshold.",
            optional = true,
            minValue = 0.,
            maxValue = 1.)
    private Double calibrationSensitivityThreshold;

    @Argument(
            fullName = MODE_LONG_NAME,
            doc = "Variant types for which to train models. Duplicate values will be ignored.",
            minElements = 1)
    public List<VariantType> variantTypes = new ArrayList<>(Arrays.asList(VariantType.SNP, VariantType.INDEL));

    private AvailableLabelsMode availableLabelsMode;

    @Override
    protected Object doWork() {

        validateArgumentsAndSetModes();

        logger.info("Starting training...");

        for (final VariantType variantType : VariantType.values()) { // enforces order in which models are trained
            if (variantTypes.contains(variantType)) {
                doModelingWorkForVariantType(variantType);
            }
        }

        logger.info(String.format("%s complete.", getClass().getSimpleName()));

        return null;
    }

    private void validateArgumentsAndSetModes() {
        IOUtils.canReadFile(inputAnnotationsFile);

        Utils.validateArg((inputUnlabeledAnnotationsFile == null) == (calibrationSensitivityThreshold == null),
                "Unlabeled annotations and calibration-sensitivity threshold must both be unspecified (for positive-only model training) " +
                        "or specified (for positive-unlabeled model training).");

        availableLabelsMode = inputUnlabeledAnnotationsFile != null && calibrationSensitivityThreshold != null
                ? AvailableLabelsMode.POSITIVE_UNLABELED
                : AvailableLabelsMode.POSITIVE_ONLY;

        if (inputUnlabeledAnnotationsFile != null) {
            IOUtils.canReadFile(inputUnlabeledAnnotationsFile);
            final List<String> annotationNames = LabeledVariantAnnotationsData.readAnnotationNames(inputAnnotationsFile);
            final List<String> unlabeledAnnotationNames = LabeledVariantAnnotationsData.readAnnotationNames(inputUnlabeledAnnotationsFile);
            Utils.validateArg(annotationNames.equals(unlabeledAnnotationNames), "Annotation names must be identical for positive and unlabeled annotations.");
        }

        switch (modelBackend) {
            case JAVA_BGMM:
                Utils.validateArg(pythonScriptFile == null,
                        "Python script should not be provided when using JAVA_BGMM backend.");
                IOUtils.canReadFile(hyperparametersJSONFile);
                logger.info("Running in JAVA_BGMM mode...");
                break;
            case PYTHON_IFOREST:
                Utils.validateArg(pythonScriptFile == null,
                        "Python script should not be provided when using PYTHON_IFOREST backend.");

                pythonScriptFile = IOUtils.writeTempResource(new Resource(ISOLATION_FOREST_PYTHON_SCRIPT, TrainVariantAnnotationsModel.class));
                if (hyperparametersJSONFile == null) {
                    hyperparametersJSONFile = IOUtils.writeTempResource(new Resource(ISOLATION_FOREST_HYPERPARAMETERS_JSON, TrainVariantAnnotationsModel.class));
                }
                IOUtils.canReadFile(hyperparametersJSONFile);
                PythonScriptExecutor.checkPythonEnvironmentForPackage("argparse");
                PythonScriptExecutor.checkPythonEnvironmentForPackage("h5py");
                PythonScriptExecutor.checkPythonEnvironmentForPackage("numpy");
                PythonScriptExecutor.checkPythonEnvironmentForPackage("sklearn");
                PythonScriptExecutor.checkPythonEnvironmentForPackage("dill");
                logger.info("Running in PYTHON_IFOREST mode...");
                break;
            case PYTHON_SCRIPT:
                IOUtils.canReadFile(pythonScriptFile);
                IOUtils.canReadFile(hyperparametersJSONFile);
                logger.info("Running in PYTHON_SCRIPT mode...");
                break;
            default:
                throw new GATKException.ShouldNeverReachHereException("Unknown model-backend mode.");
        }
    }

    /**
     * TODO
     */
    private void doModelingWorkForVariantType(final VariantType variantType) {
        // positive model
        final List<String> annotationNames = LabeledVariantAnnotationsData.readAnnotationNames(inputAnnotationsFile);
        final double[][] annotations = LabeledVariantAnnotationsData.readAnnotations(inputAnnotationsFile);

        final List<Boolean> isTraining = LabeledVariantAnnotationsData.readLabel(inputAnnotationsFile, LabeledVariantAnnotationsData.TRAINING_LABEL);
        final List<Boolean> isCalibration = LabeledVariantAnnotationsData.readLabel(inputAnnotationsFile, LabeledVariantAnnotationsData.CALIBRATION_LABEL);
        final List<Boolean> isSNP = LabeledVariantAnnotationsData.readLabel(inputAnnotationsFile, LabeledVariantAnnotationsData.SNP_LABEL);
        final List<Boolean> isVariantType = variantType == VariantType.SNP ? isSNP : isSNP.stream().map(x -> !x).collect(Collectors.toList());

        final List<Boolean> isTrainingAndVariantType = Streams.zip(isTraining.stream(), isVariantType.stream(), (a, b) -> a && b).collect(Collectors.toList());
        final int numTrainingAndVariantType = numPassingFilter(isTrainingAndVariantType);

        final String variantTypeString = variantType.toString();
        final String outputPrefixTag = '.' + variantType.toString().toLowerCase();

        if (numTrainingAndVariantType > 0) {
            logger.info(String.format("Training %s model with %d training sites x %d annotations %s...",
                    variantTypeString, numTrainingAndVariantType, annotationNames.size(), annotationNames));
            final File labeledTrainingAndVariantTypeAnnotationsFile = LabeledVariantAnnotationsData.subsetAnnotationsToTemporaryFile(annotationNames, annotations, isTrainingAndVariantType);
            trainAndSerializeModel(labeledTrainingAndVariantTypeAnnotationsFile, outputPrefixTag);
            logger.info(String.format("%s model trained and serialized with output prefix \"%s\".", variantTypeString, outputPrefix + outputPrefixTag));

            if (modelBackend == VariantAnnotationsModelBackend.JAVA_BGMM) {
                BGMMVariantAnnotationsScorer.preprocessAnnotationsWithBGMMAndWriteHDF5(
                        annotationNames, outputPrefix + outputPrefixTag, labeledTrainingAndVariantTypeAnnotationsFile, logger);
            }

            logger.info(String.format("Scoring %d %s training sites...", numTrainingAndVariantType, variantTypeString));
            final File labeledTrainingAndVariantTypeScoresFile = score(labeledTrainingAndVariantTypeAnnotationsFile, outputPrefixTag, TRAINING_SCORES_HDF5_SUFFIX);
            logger.info(String.format("%s training scores written to %s.", variantTypeString, labeledTrainingAndVariantTypeScoresFile.getAbsolutePath()));

            final List<Boolean> isLabeledCalibrationAndVariantType = Streams.zip(isCalibration.stream(), isVariantType.stream(), (a, b) -> a && b).collect(Collectors.toList());
            final int numLabeledCalibrationAndVariantType = numPassingFilter(isLabeledCalibrationAndVariantType);
            if (numLabeledCalibrationAndVariantType > 0) {
                logger.info(String.format("Scoring %d %s calibration sites...", numLabeledCalibrationAndVariantType, variantTypeString));
                final File labeledCalibrationAndVariantTypeAnnotationsFile = LabeledVariantAnnotationsData.subsetAnnotationsToTemporaryFile(annotationNames, annotations, isLabeledCalibrationAndVariantType);
                final File labeledCalibrationAndVariantTypeScoresFile = score(labeledCalibrationAndVariantTypeAnnotationsFile, outputPrefixTag, CALIBRATION_SCORES_HDF5_SUFFIX);
                logger.info(String.format("%s calibration scores written to %s.", variantTypeString, labeledCalibrationAndVariantTypeScoresFile.getAbsolutePath()));
            } else {
                logger.warn(String.format("No %s calibration sites were available.", variantTypeString));
            }

            // negative model
            if (availableLabelsMode == AvailableLabelsMode.POSITIVE_UNLABELED) {
                final double[][] unlabeledAnnotations = LabeledVariantAnnotationsData.readAnnotations(inputUnlabeledAnnotationsFile);
                final List<Boolean> unlabeledIsSNP = LabeledVariantAnnotationsData.readLabel(inputUnlabeledAnnotationsFile, "snp");
                final List<Boolean> isUnlabeledVariantType = variantType == VariantType.SNP ? unlabeledIsSNP : unlabeledIsSNP.stream().map(x -> !x).collect(Collectors.toList());

                final int numUnlabeledVariantType = numPassingFilter(isUnlabeledVariantType);

                if (numUnlabeledVariantType > 0) {
                    final File labeledCalibrationAndVariantTypeScoresFile = new File(outputPrefix + outputPrefixTag + CALIBRATION_SCORES_HDF5_SUFFIX);
                    final double[] labeledCalibrationAndVariantTypeScores = VariantAnnotationsScorer.readScores(labeledCalibrationAndVariantTypeScoresFile);
                    final double scoreThreshold = calibrationSensitivityThreshold == 1. // Percentile requires quantile > 0, so we treat this as a special case
                            ? Doubles.min(labeledCalibrationAndVariantTypeScores)
                            : new Percentile(100. * (1. - calibrationSensitivityThreshold)).evaluate(labeledCalibrationAndVariantTypeScores);
                    logger.info(String.format("Using %s score threshold of %.4f corresponding to specified calibration-sensitivity threshold of %.4f ...",
                            variantTypeString, scoreThreshold, calibrationSensitivityThreshold));

                    final double[] labeledTrainingAndVariantTypeScores = VariantAnnotationsScorer.readScores(labeledTrainingAndVariantTypeScoresFile);
                    final List<Boolean> isNegativeTrainingFromLabeledTrainingAndVariantType = Arrays.stream(labeledTrainingAndVariantTypeScores).boxed().map(s -> s < scoreThreshold).collect(Collectors.toList());
                    final int numNegativeTrainingFromLabeledTrainingAndVariantType = numPassingFilter(isNegativeTrainingFromLabeledTrainingAndVariantType);
                    logger.info(String.format("Selected %d labeled %s sites below score threshold of %.4f for negative-model training...",
                            numNegativeTrainingFromLabeledTrainingAndVariantType, variantTypeString, scoreThreshold));

                    logger.info(String.format("Scoring %d unlabeled %s sites...", numUnlabeledVariantType, variantTypeString));
                    final File unlabeledVariantTypeAnnotationsFile = LabeledVariantAnnotationsData.subsetAnnotationsToTemporaryFile(annotationNames, unlabeledAnnotations, isUnlabeledVariantType);
                    final File unlabeledVariantTypeScoresFile = score(unlabeledVariantTypeAnnotationsFile, outputPrefixTag, UNLABELED_SCORES_HDF5_SUFFIX);
                    final double[] unlabeledVariantTypeScores = VariantAnnotationsScorer.readScores(unlabeledVariantTypeScoresFile);
                    final List<Boolean> isNegativeTrainingFromUnlabeledVariantType = Arrays.stream(unlabeledVariantTypeScores).boxed().map(s -> s < scoreThreshold).collect(Collectors.toList()); // length matches unlabeledAnnotationsFile
                    final int numNegativeTrainingFromUnlabeledVariantType = numPassingFilter(isNegativeTrainingFromUnlabeledVariantType);
                    logger.info(String.format("Selected %d unlabeled %s sites below score threshold of %.4f for negative-model training...",
                            numNegativeTrainingFromUnlabeledVariantType, variantTypeString, scoreThreshold));

                    final double[][] negativeTrainingAndVariantTypeAnnotations = concatenateLabeledAndUnlabeledNegativeTrainingData(
                            annotationNames, annotations, unlabeledAnnotations, isNegativeTrainingFromLabeledTrainingAndVariantType, isNegativeTrainingFromUnlabeledVariantType);
                    final int numNegativeTrainingAndVariantType = negativeTrainingAndVariantTypeAnnotations.length;
                    final List<Boolean> isNegativeTrainingAndVariantType = Collections.nCopies(numNegativeTrainingAndVariantType, true);

                    logger.info(String.format("Training %s negative model with %d negative-training sites x %d annotations %s...",
                            variantTypeString, numNegativeTrainingAndVariantType, annotationNames.size(), annotationNames));
                    final File negativeTrainingAnnotationsFile = LabeledVariantAnnotationsData.subsetAnnotationsToTemporaryFile(
                            annotationNames, negativeTrainingAndVariantTypeAnnotations, isNegativeTrainingAndVariantType);
                    trainAndSerializeModel(negativeTrainingAnnotationsFile, outputPrefixTag + NEGATIVE_TAG);
                    logger.info(String.format("%s negative model trained and serialized with output prefix \"%s\".", variantTypeString, outputPrefix + outputPrefixTag + NEGATIVE_TAG));

                    if (numLabeledCalibrationAndVariantType > 0) {
                        logger.info(String.format("Re-scoring %d %s calibration sites...", numLabeledCalibrationAndVariantType, variantTypeString));
                        final File labeledCalibrationAnnotationsFile = LabeledVariantAnnotationsData.subsetAnnotationsToTemporaryFile(annotationNames, annotations, isLabeledCalibrationAndVariantType);
                        final File labeledCalibrationScoresFile = positiveNegativeScore(labeledCalibrationAnnotationsFile, outputPrefixTag, CALIBRATION_SCORES_HDF5_SUFFIX);
                        logger.info(String.format("Calibration scores written to %s.", labeledCalibrationScoresFile.getAbsolutePath()));
                    }
                } else {
                    throw new UserException.BadInput(String.format("Attempted to train %s negative model, " +
                            "but no suitable sites were found in the provided annotations.", variantTypeString));
                }
            }
        } else {
            throw new UserException.BadInput(String.format("Attempted to train %s model, " +
                    "but no suitable training sites were found in the provided annotations.", variantTypeString));
        }
    }

    private static int numPassingFilter(List<Boolean> isPassing) {
        return isPassing.stream().mapToInt(x -> x ? 1 : 0).sum();
    }

    private void trainAndSerializeModel(final File trainingAnnotationsFile,
                                        final String outputPrefixTag) {
        readAndValidateTrainingAnnotations(trainingAnnotationsFile, outputPrefixTag);
        final VariantAnnotationsModel model;
        switch (modelBackend) {
            case JAVA_BGMM:
                model = new BGMMVariantAnnotationsModel(hyperparametersJSONFile);
                break;
            case PYTHON_IFOREST:
                model = new PythonSklearnVariantAnnotationsModel(pythonScriptFile, hyperparametersJSONFile);
                break;
            case PYTHON_SCRIPT:
                model = new PythonSklearnVariantAnnotationsModel(pythonScriptFile, hyperparametersJSONFile);
                break;
            default:
                throw new GATKException.ShouldNeverReachHereException("Unknown model mode.");
        }
        model.trainAndSerialize(trainingAnnotationsFile, outputPrefix + outputPrefixTag);
    }

    /**
     * When training models on data that has been subset to a given variant type,
     * we FAIL if any annotation is completely missing and WARN if any annotation has zero variance.
     */
    private void readAndValidateTrainingAnnotations(final File trainingAnnotationsFile,
                                                    final String outputPrefixTag) {
        final List<String> annotationNames = LabeledVariantAnnotationsData.readAnnotationNames(trainingAnnotationsFile);
        final double[][] annotations = LabeledVariantAnnotationsData.readAnnotations(trainingAnnotationsFile);

        // these checks are redundant, but we err on the side of robustness
        final int numAnnotationNames = annotationNames.size();
        final int numData = annotations.length;
        Utils.validateArg(numAnnotationNames > 0, "Number of annotation names must be positive.");
        Utils.validateArg(numData > 0, "Number of data points must be positive.");
        final int numFeatures = annotations[0].length;
        Utils.validateArg(numAnnotationNames == numFeatures,
                "Number of annotation names must match the number of features in the annotation data.");

        final List<String> completelyMissingAnnotationNames = new ArrayList<>(numFeatures);
        IntStream.range(0, numFeatures).forEach(
                i -> {
                    if (new Variance().evaluate(IntStream.range(0, numData).mapToDouble(n -> annotations[n][i]).toArray()) == 0.) {
                        logger.warn(String.format("All values of the annotation %s are identical in the training data for the %s model.",
                                annotationNames.get(i), outputPrefix + outputPrefixTag));
                    }
                    if (IntStream.range(0, numData).boxed().map(n -> annotations[n][i]).allMatch(x -> Double.isNaN(x))) {
                        completelyMissingAnnotationNames.add(annotationNames.get(i));
                    }
                }
        );

        if (!completelyMissingAnnotationNames.isEmpty()) {
            throw new UserException.BadInput(
                    String.format("All values of the following annotations are missing in the training data for the %s model: %s. " +
                                    "Consider repeating the extraction step with this annotation dropped. " +
                                    "If this is a negative model and the amount of negative training data is small, " +
                                    "perhaps also consider lowering the value of the %s argument so that more " +
                                    "training data is considered, which may ultimately admit data with non-missing values for the annotation " +
                                    "(although note that this will also have implications for the resulting model fit); " +
                                    "alternatively, consider excluding the %s and %s arguments and running positive-only modeling.",
                            outputPrefix + outputPrefixTag, completelyMissingAnnotationNames,
                            CALIBRATION_SENSITIVITY_THRESHOLD_LONG_NAME, UNLABELED_ANNOTATIONS_HDF5_LONG_NAME, CALIBRATION_SENSITIVITY_THRESHOLD_LONG_NAME));
        }
    }

    private File score(final File annotationsFile,
                       final String outputPrefixTag,
                       final String outputSuffix) {
        final VariantAnnotationsScorer scorer;
        switch (modelBackend) {
            case JAVA_BGMM:
                scorer = BGMMVariantAnnotationsScorer.deserialize(new File(outputPrefix + outputPrefixTag + BGMMVariantAnnotationsScorer.BGMM_SCORER_SER_SUFFIX));
                break;
            case PYTHON_IFOREST:
            case PYTHON_SCRIPT:
                scorer = new PythonSklearnVariantAnnotationsScorer(pythonScriptFile, new File(outputPrefix + outputPrefixTag + PythonSklearnVariantAnnotationsScorer.PYTHON_SCORER_PKL_SUFFIX));
                break;

            default:
                throw new GATKException.ShouldNeverReachHereException("Unknown model mode.");
        }
        final File outputScoresFile = new File(outputPrefix + outputPrefixTag + outputSuffix);
        scorer.score(annotationsFile, outputScoresFile);
        return outputScoresFile;
    }

    private File positiveNegativeScore(final File annotationsFile,
                                       final String outputPrefixTag,
                                       final String outputSuffix) {
        final VariantAnnotationsScorer scorer;
        switch (modelBackend) {
            case JAVA_BGMM:
                scorer = VariantAnnotationsScorer.combinePositiveAndNegativeScorer(
                        BGMMVariantAnnotationsScorer.deserialize(new File(outputPrefix + outputPrefixTag + BGMMVariantAnnotationsScorer.BGMM_SCORER_SER_SUFFIX)),
                        BGMMVariantAnnotationsScorer.deserialize(new File(outputPrefix + outputPrefixTag + NEGATIVE_TAG + BGMMVariantAnnotationsScorer.BGMM_SCORER_SER_SUFFIX)));
                break;
            case PYTHON_IFOREST:
            case PYTHON_SCRIPT:
                scorer = VariantAnnotationsScorer.combinePositiveAndNegativeScorer(
                        new PythonSklearnVariantAnnotationsScorer(pythonScriptFile, new File(outputPrefix + outputPrefixTag + PythonSklearnVariantAnnotationsScorer.PYTHON_SCORER_PKL_SUFFIX)),
                        new PythonSklearnVariantAnnotationsScorer(pythonScriptFile, new File(outputPrefix + outputPrefixTag + NEGATIVE_TAG + PythonSklearnVariantAnnotationsScorer.PYTHON_SCORER_PKL_SUFFIX)));
                break;
            default:
                throw new GATKException.ShouldNeverReachHereException("Unknown model mode.");
        }
        final File outputScoresFile = new File(outputPrefix + outputPrefixTag + outputSuffix);
        scorer.score(annotationsFile, outputScoresFile);
        return outputScoresFile;
    }

    private static double[][] concatenateLabeledAndUnlabeledNegativeTrainingData(final List<String> annotationNames,
                                                                                 final double[][] annotations,
                                                                                 final double[][] unlabeledAnnotations,
                                                                                 final List<Boolean> isNegativeTrainingFromLabeledTrainingAndVariantType,
                                                                                 final List<Boolean> isNegativeTrainingFromUnlabeledVariantType) {
        final File negativeTrainingFromLabeledTrainingAndVariantTypeAnnotationsFile =
                LabeledVariantAnnotationsData.subsetAnnotationsToTemporaryFile(annotationNames, annotations, isNegativeTrainingFromLabeledTrainingAndVariantType);
        final double[][] negativeTrainingFromLabeledTrainingAndVariantTypeAnnotations = LabeledVariantAnnotationsData.readAnnotations(negativeTrainingFromLabeledTrainingAndVariantTypeAnnotationsFile);

        final File negativeTrainingFromUnlabeledVariantTypeAnnotationsFile =
                LabeledVariantAnnotationsData.subsetAnnotationsToTemporaryFile(annotationNames, unlabeledAnnotations, isNegativeTrainingFromUnlabeledVariantType);
        final double[][] negativeTrainingFromUnlabeledVariantTypeAnnotations = LabeledVariantAnnotationsData.readAnnotations(negativeTrainingFromUnlabeledVariantTypeAnnotationsFile);

        return Streams.concat(
                Arrays.stream(negativeTrainingFromLabeledTrainingAndVariantTypeAnnotations),
                Arrays.stream(negativeTrainingFromUnlabeledVariantTypeAnnotations)).toArray(double[][]::new);
    }
}