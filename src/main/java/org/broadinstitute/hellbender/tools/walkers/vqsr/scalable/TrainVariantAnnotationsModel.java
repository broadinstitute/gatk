package org.broadinstitute.hellbender.tools.walkers.vqsr.scalable;

import com.google.common.collect.Streams;
import com.google.common.primitives.Doubles;
import org.apache.commons.math3.stat.descriptive.rank.Percentile;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hdf5.HDF5File;
import org.broadinstitute.hdf5.HDF5LibException;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.utils.HDF5Utils;
import org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.data.LabeledVariantAnnotationsData;
import org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.data.VariantType;
import org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.modeling.BGMMVariantAnnotationsModel;
import org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.modeling.BGMMVariantAnnotationsScorer;
import org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.modeling.PythonSklearnVariantAnnotationsModel;
import org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.modeling.PythonSklearnVariantAnnotationsScorer;
import org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.modeling.VariantAnnotationsModel;
import org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.modeling.VariantAnnotationsScorer;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.python.PythonScriptExecutor;
import picard.cmdline.programgroups.VariantFilteringProgramGroup;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

/**
 * TODO
 */
@CommandLineProgramProperties(
        // TODO
        summary = "",
        oneLineSummary = "",
        programGroup = VariantFilteringProgramGroup.class
)
@DocumentedFeature
public final class TrainVariantAnnotationsModel extends CommandLineProgram {

    public static final String MODE_LONG_NAME = "mode";
    public static final String ANNOTATIONS_HDF5_LONG_NAME = "annotations-hdf5";
    public static final String UNLABELED_ANNOTATIONS_HDF5_LONG_NAME = "unlabeled-annotations-hdf5";
    public static final String PYTHON_SCRIPT_LONG_NAME = "python-script";
    public static final String HYPERPARAMETERS_JSON_LONG_NAME = "hyperparameters-json";
    public static final String CALIBRATION_SENSITIVITY_THRESHOLD_LONG_NAME = "calibration-sensitivity-threshold";

    enum ModelBackendMode {
        PYTHON, BGMM    // TODO put IsolationForest script into resources and use as a default PYTHON backend
    }

    enum AvailableLabelsMode {
        POSITIVE_ONLY, POSITIVE_UNLABELED
    }

    private static final String TRAINING_SCORES_HDF5_SUFFIX = ".trainingScores.hdf5";
    private static final String CALIBRATION_SCORES_HDF5_SUFFIX = ".calibrationScores.hdf5";
    private static final String UNLABELED_SCORES_HDF5_SUFFIX = ".unlabeledScores.hdf5";

    @Argument(
            fullName = ANNOTATIONS_HDF5_LONG_NAME,
            doc = "HDF5 file containing annotations extracted with ExtractVariantAnnotations.")
    private File inputAnnotationsFile;

    @Argument(
            fullName = UNLABELED_ANNOTATIONS_HDF5_LONG_NAME,
            optional = true,
            doc = "HDF5 file containing annotations extracted with ExtractVariantAnnotations.")
    private File inputUnlabeledAnnotationsFile;

    @Argument(
            fullName = PYTHON_SCRIPT_LONG_NAME,
            optional = true)
    private File pythonScriptFile;

    @Argument(
            fullName = HYPERPARAMETERS_JSON_LONG_NAME,
            doc = "JSON file containing hyperparameters.")
    private File hyperparametersJSONFile;

    @Argument(
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "Output prefix.")
    private String outputPrefix;

    @Argument(
            fullName = CALIBRATION_SENSITIVITY_THRESHOLD_LONG_NAME,
            doc = "Calibration-sensitivity threshold that determines which sites will be used for training the negative model. " +
                    "Increasing this will decrease the corresponding positive-model score threshold; sites with scores below this score " +
                    "threshold will be used for training the negative model. Thus, this parameter should typically be chosen to " +
                    "be close to 1, so that sites that score highly according to the positive model will not be used to train the negative model.",
            optional = true,
            minValue = 0.,
            maxValue = 1.)
    private Double calibrationSensitivityThreshold;

    @Argument(
            fullName = MODE_LONG_NAME,
            doc = "Variant types for which to train models. Duplicate values will be ignored.",
            minElements = 1,
            optional = true)
    public List<VariantType> variantTypes = new ArrayList<>(Arrays.asList(VariantType.SNP, VariantType.INDEL));

    private ModelBackendMode modelBackendMode;
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
        IOUtils.canReadFile(hyperparametersJSONFile);

        Utils.validateArg((inputUnlabeledAnnotationsFile == null) == (calibrationSensitivityThreshold == null),
                "Unlabeled annotations and calibration-sensitivity threshold must both be unspecified (for positive-only model training) " +
                        "or specified (for positive-negative model training).");

        availableLabelsMode = inputUnlabeledAnnotationsFile != null && calibrationSensitivityThreshold != null
                ? AvailableLabelsMode.POSITIVE_UNLABELED
                : AvailableLabelsMode.POSITIVE_ONLY;

        if (inputUnlabeledAnnotationsFile != null) {
            IOUtils.canReadFile(inputUnlabeledAnnotationsFile);
            final List<String> annotationNames = LabeledVariantAnnotationsData.readAnnotationNames(inputAnnotationsFile);
            final List<String> unlabeledAnnotationNames = LabeledVariantAnnotationsData.readAnnotationNames(inputUnlabeledAnnotationsFile);
            Utils.validateArg(annotationNames.equals(unlabeledAnnotationNames), "Annotation names must be identical for positive and unlabeled annotations.");
        }

        // TODO test output and fail early
        // TODO FAIL if annotations that are all NaN?
        // TODO WARN if annotations that have zero variance?

        if (pythonScriptFile != null) {
            logger.info("Python script was provided, running in PYTHON mode...");
            modelBackendMode = ModelBackendMode.PYTHON;

            IOUtils.canReadFile(pythonScriptFile);
            PythonScriptExecutor.checkPythonEnvironmentForPackage("argparse");
            PythonScriptExecutor.checkPythonEnvironmentForPackage("h5py");
            PythonScriptExecutor.checkPythonEnvironmentForPackage("numpy");
            PythonScriptExecutor.checkPythonEnvironmentForPackage("sklearn");
            PythonScriptExecutor.checkPythonEnvironmentForPackage("dill");
        } else {
            logger.info("Python script was not provided, running in BGMM mode...");
            modelBackendMode = ModelBackendMode.BGMM;
        }
    }

    private void doModelingWorkForVariantType(final VariantType variantType) {
        // positive model
        final List<String> annotationNames = LabeledVariantAnnotationsData.readAnnotationNames(inputAnnotationsFile);
        final double[][] annotations = LabeledVariantAnnotationsData.readAnnotations(inputAnnotationsFile);

        final List<Boolean> isTraining = LabeledVariantAnnotationsData.readLabel(inputAnnotationsFile, LabeledVariantAnnotationsData.TRAINING_LABEL);
        final List<Boolean> isCalibration = LabeledVariantAnnotationsData.readLabel(inputAnnotationsFile, LabeledVariantAnnotationsData.CALIBRATION_LABEL);
        final List<Boolean> isSNP = LabeledVariantAnnotationsData.readLabel(inputAnnotationsFile, "snp");
        final List<Boolean> isVariantType = variantType == VariantType.SNP ? isSNP : isSNP.stream().map(x -> !x).collect(Collectors.toList());

        final List<Boolean> isTrainingAndVariantType = Streams.zip(isTraining.stream(), isVariantType.stream(), (a, b) -> a && b).collect(Collectors.toList());
        final int numTrainingAndVariantType = numPassingFilter(isTrainingAndVariantType);

        final String logMessageTag = variantType.toString();
        final String outputPrefixTag = '.' + variantType.toString().toLowerCase();

        if (numTrainingAndVariantType > 0) {
            logger.info(String.format("Training %s model with %d training sites x %d annotations %s...",
                    logMessageTag, numTrainingAndVariantType, annotationNames.size(), annotationNames));
            final File labeledTrainingAndVariantTypeAnnotationsFile = LabeledVariantAnnotationsData.subsetAnnotationsToTemporaryFile(annotationNames, annotations, isTrainingAndVariantType);
            trainAndSerializeModel(labeledTrainingAndVariantTypeAnnotationsFile, outputPrefixTag);
            logger.info(String.format("%s model trained and serialized with output prefix \"%s\".", logMessageTag, outputPrefix + outputPrefixTag));

            if (modelBackendMode == ModelBackendMode.BGMM) {
                preprocessAnnotationsWithBGMMAndWriteHDF5(annotationNames, outputPrefixTag, labeledTrainingAndVariantTypeAnnotationsFile);
            }

            logger.info(String.format("Scoring %d %s training sites...", numTrainingAndVariantType, logMessageTag));
            final File labeledTrainingAndVariantTypeScoresFile = score(labeledTrainingAndVariantTypeAnnotationsFile, outputPrefixTag, TRAINING_SCORES_HDF5_SUFFIX);
            logger.info(String.format("%s training scores written to %s.", logMessageTag, labeledTrainingAndVariantTypeScoresFile.getAbsolutePath()));

            final List<Boolean> isLabeledCalibrationAndVariantType = Streams.zip(isCalibration.stream(), isVariantType.stream(), (a, b) -> a && b).collect(Collectors.toList());
            final int numLabeledCalibrationAndVariantType = numPassingFilter(isLabeledCalibrationAndVariantType);
            if (numLabeledCalibrationAndVariantType > 0) {
                logger.info(String.format("Scoring %d %s calibration sites...", numLabeledCalibrationAndVariantType, logMessageTag));
                final File labeledCalibrationAndVariantTypeAnnotationsFile = LabeledVariantAnnotationsData.subsetAnnotationsToTemporaryFile(annotationNames, annotations, isLabeledCalibrationAndVariantType);
                final File labeledCalibrationAndVariantTypeScoresFile = score(labeledCalibrationAndVariantTypeAnnotationsFile, outputPrefixTag, CALIBRATION_SCORES_HDF5_SUFFIX);
                logger.info(String.format("%s calibration scores written to %s.", logMessageTag, labeledCalibrationAndVariantTypeScoresFile.getAbsolutePath()));
            }

            // negative model
            if (availableLabelsMode == AvailableLabelsMode.POSITIVE_UNLABELED) {
                final List<String> unlabeledAnnotationNames = LabeledVariantAnnotationsData.readAnnotationNames(inputUnlabeledAnnotationsFile);
                final double[][] unlabeledAnnotations = LabeledVariantAnnotationsData.readAnnotations(inputUnlabeledAnnotationsFile);
                final List<Boolean> unlabeledIsSNP = LabeledVariantAnnotationsData.readLabel(inputUnlabeledAnnotationsFile, "snp");
                final List<Boolean> isUnlabeledVariantType = variantType == VariantType.SNP ? unlabeledIsSNP : unlabeledIsSNP.stream().map(x -> !x).collect(Collectors.toList());

                final int numUnlabeledVariantType = numPassingFilter(isUnlabeledVariantType);

                if (numUnlabeledVariantType > 0) {
                    final File labeledCalibrationAndVariantTypeScoresFile = new File(outputPrefix + outputPrefixTag + CALIBRATION_SCORES_HDF5_SUFFIX); // produced by doModelingAndScoringWork TODO output a copy of these?
                    final double[] labeledCalibrationAndVariantTypeScores = VariantAnnotationsScorer.readScores(labeledCalibrationAndVariantTypeScoresFile);
                    final double scoreThreshold = calibrationSensitivityThreshold == 1. // Percentile requires quantile > 0, so we treat this as a special case
                            ? Doubles.min(labeledCalibrationAndVariantTypeScores)
                            : new Percentile(100. * (1. - calibrationSensitivityThreshold)).evaluate(labeledCalibrationAndVariantTypeScores);
                    logger.info(String.format("Using %s score threshold of %.4f corresponding to specified calibration-sensitivity threshold of %.4f ...",
                            logMessageTag, scoreThreshold, calibrationSensitivityThreshold));

                    final double[] labeledTrainingAndVariantTypeScores = VariantAnnotationsScorer.readScores(labeledTrainingAndVariantTypeScoresFile);
                    final List<Boolean> isNegativeTrainingFromLabeledTrainingAndVariantType = Arrays.stream(labeledTrainingAndVariantTypeScores).boxed().map(s -> s < scoreThreshold).collect(Collectors.toList());
                    final int numNegativeTrainingFromLabeledTrainingAndVariantType = numPassingFilter(isNegativeTrainingFromLabeledTrainingAndVariantType);
                    logger.info(String.format("Selected %d labeled %s sites below score threshold of %.4f for negative-model training...",
                            numNegativeTrainingFromLabeledTrainingAndVariantType, logMessageTag, scoreThreshold));

                    logger.info(String.format("Scoring %d unlabeled %s sites...", numUnlabeledVariantType, logMessageTag));
                    final File unlabeledVariantTypeAnnotationsFile = LabeledVariantAnnotationsData.subsetAnnotationsToTemporaryFile(annotationNames, unlabeledAnnotations, isUnlabeledVariantType);
                    final File unlabeledVariantTypeScoresFile = score(unlabeledVariantTypeAnnotationsFile, outputPrefixTag, UNLABELED_SCORES_HDF5_SUFFIX);
                    final double[] unlabeledVariantTypeScores = VariantAnnotationsScorer.readScores(unlabeledVariantTypeScoresFile);
                    final List<Boolean> isNegativeTrainingFromUnlabeledVariantType = Arrays.stream(unlabeledVariantTypeScores).boxed().map(s -> s < scoreThreshold).collect(Collectors.toList()); // length matches unlabeledAnnotationsFile
                    final int numNegativeTrainingFromUnlabeledVariantType = numPassingFilter(isNegativeTrainingFromUnlabeledVariantType);
                    logger.info(String.format("Selected %d unlabeled %s sites below score threshold of %.4f for negative-model training...",
                            numNegativeTrainingFromUnlabeledVariantType, logMessageTag, scoreThreshold));

                    // combine labeled and unlabeled negative training data
                    final File negativeTrainingFromLabeledTrainingAndVariantTypeAnnotationsFile =
                            LabeledVariantAnnotationsData.subsetAnnotationsToTemporaryFile(annotationNames, annotations, isNegativeTrainingFromLabeledTrainingAndVariantType);
                    final double[][] negativeTrainingFromLabeledTrainingAndVariantTypeAnnotations = LabeledVariantAnnotationsData.readAnnotations(negativeTrainingFromLabeledTrainingAndVariantTypeAnnotationsFile);

                    final File negativeTrainingFromUnlabeledVariantTypeAnnotationsFile =
                            LabeledVariantAnnotationsData.subsetAnnotationsToTemporaryFile(annotationNames, unlabeledAnnotations, isNegativeTrainingFromUnlabeledVariantType);
                    final double[][] negativeTrainingFromUnlabeledVariantTypeAnnotations = LabeledVariantAnnotationsData.readAnnotations(negativeTrainingFromUnlabeledVariantTypeAnnotationsFile);

                    final double[][] negativeTrainingAndVariantTypeAnnotations = Streams.concat(
                            Arrays.stream(negativeTrainingFromLabeledTrainingAndVariantTypeAnnotations),
                            Arrays.stream(negativeTrainingFromUnlabeledVariantTypeAnnotations)).toArray(double[][]::new);
                    final int numNegativeTrainingAndVariantType = negativeTrainingAndVariantTypeAnnotations.length;
                    final List<Boolean> isNegativeTrainingAndVariantType = Collections.nCopies(numNegativeTrainingAndVariantType, true);

                    logger.info(String.format("Training %s negative model with %d negative-training sites x %d annotations %s...",
                            logMessageTag, numNegativeTrainingAndVariantType, annotationNames.size(), annotationNames));
                    final File negativeTrainingAnnotationsFile = LabeledVariantAnnotationsData.subsetAnnotationsToTemporaryFile(
                            annotationNames, negativeTrainingAndVariantTypeAnnotations, isNegativeTrainingAndVariantType);
                    trainAndSerializeModel(negativeTrainingAnnotationsFile, outputPrefixTag + ".negative");
                    logger.info(String.format("%s negative model trained and serialized with output prefix \"%s\".", logMessageTag, outputPrefix + outputPrefixTag + ".negative"));

                    if (numLabeledCalibrationAndVariantType > 0) {
                        logger.info(String.format("Re-scoring %d %s calibration sites...", numLabeledCalibrationAndVariantType, logMessageTag));
                        final File labeledCalibrationAnnotationsFile = LabeledVariantAnnotationsData.subsetAnnotationsToTemporaryFile(annotationNames, annotations, isLabeledCalibrationAndVariantType);
                        final File labeledCalibrationScoresFile = positiveNegativeScore(labeledCalibrationAnnotationsFile, outputPrefixTag, CALIBRATION_SCORES_HDF5_SUFFIX);
                        logger.info(String.format("Calibration scores written to %s.", labeledCalibrationScoresFile.getAbsolutePath()));
                    }
                } else {
                    throw new UserException.BadInput(String.format("Attempted to train %s negative model, " +
                            "but no suitable sites were found in the provided annotations.", logMessageTag));
                }
            }
        } else {
            throw new UserException.BadInput(String.format("Attempted to train %s model, " +
                    "but no suitable training sites were found in the provided annotations.", logMessageTag));
        }
    }

    private static int numPassingFilter(List<Boolean> isPassing) {
        return isPassing.stream().mapToInt(x -> x ? 1 : 0).sum();
    }

    // TODO clean this up
    private void preprocessAnnotationsWithBGMMAndWriteHDF5(final List<String> annotationNames,
                                                           final String outputPrefixTag,
                                                           final File labeledTrainingAndVariantTypeAnnotationsFile) {
        final double[][] rawAnnotations = LabeledVariantAnnotationsData.readAnnotations(labeledTrainingAndVariantTypeAnnotationsFile);
        final BGMMVariantAnnotationsScorer scorer = BGMMVariantAnnotationsScorer.deserialize(new File(outputPrefix + outputPrefixTag + BGMMVariantAnnotationsModel.BGMM_SCORER_SER_SUFFIX));
        final double[][] preprocessedAnnotations = scorer.preprocess(rawAnnotations);
        final File outputPreprocessedAnnotationsFile = new File(outputPrefix + outputPrefixTag + ".annot.pre.hdf5");
        try (final HDF5File hdf5File = new HDF5File(outputPreprocessedAnnotationsFile, HDF5File.OpenMode.CREATE)) {
            IOUtils.canReadFile(hdf5File.getFile());
            hdf5File.makeStringArray("/data/annotation_names", annotationNames.toArray(new String[0]));
            HDF5Utils.writeChunkedDoubleMatrix(hdf5File, "/data/annotations", preprocessedAnnotations, HDF5Utils.MAX_NUMBER_OF_VALUES_PER_HDF5_MATRIX / 16);
        } catch (final HDF5LibException exception) {
            throw new GATKException(String.format("Exception encountered during writing of preprocessed annotations (%s). Output file at %s may be in a bad state.",
                    exception, outputPreprocessedAnnotationsFile.getAbsolutePath()));
        }
        logger.info(String.format("Preprocessed annotations written to %s.", outputPreprocessedAnnotationsFile.getAbsolutePath()));
    }

    private void trainAndSerializeModel(final File trainingAnnotationsFile,
                                        final String outputPrefixTag) {
        final VariantAnnotationsModel model;
        switch (modelBackendMode) {
            case PYTHON:
                model = new PythonSklearnVariantAnnotationsModel(pythonScriptFile, hyperparametersJSONFile);
                break;
            case BGMM:
                model = new BGMMVariantAnnotationsModel(hyperparametersJSONFile);
                break;
            default:
                throw new GATKException.ShouldNeverReachHereException("Unknown model mode.");
        }
        model.trainAndSerialize(trainingAnnotationsFile, outputPrefix + outputPrefixTag);
    }

    private File score(final File annotationsFile,
                       final String outputPrefixTag,
                       final String outputSuffix) {
        final VariantAnnotationsScorer scorer;
        switch (modelBackendMode) {
            case PYTHON:
                scorer = new PythonSklearnVariantAnnotationsScorer(pythonScriptFile, new File(outputPrefix + outputPrefixTag + ".scorer.pkl"));
                break;
            case BGMM:
                scorer = BGMMVariantAnnotationsScorer.deserialize(new File(outputPrefix + outputPrefixTag + BGMMVariantAnnotationsModel.BGMM_SCORER_SER_SUFFIX));
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
        switch (modelBackendMode) {
            case PYTHON:
                scorer = VariantAnnotationsScorer.combinePositiveAndNegativeScorer(
                        new PythonSklearnVariantAnnotationsScorer(pythonScriptFile, new File(outputPrefix + outputPrefixTag + ".scorer.pkl")),
                        new PythonSklearnVariantAnnotationsScorer(pythonScriptFile, new File(outputPrefix + outputPrefixTag + ".negative.scorer.pkl")));
                break;
            case BGMM:
                scorer = VariantAnnotationsScorer.combinePositiveAndNegativeScorer(
                        BGMMVariantAnnotationsScorer.deserialize(new File(outputPrefix + outputPrefixTag + BGMMVariantAnnotationsModel.BGMM_SCORER_SER_SUFFIX)),
                        BGMMVariantAnnotationsScorer.deserialize(new File(outputPrefix + outputPrefixTag + ".negative" + BGMMVariantAnnotationsModel.BGMM_SCORER_SER_SUFFIX)));
                break;
            default:
                throw new GATKException.ShouldNeverReachHereException("Unknown model mode.");
        }
        final File outputScoresFile = new File(outputPrefix + outputPrefixTag + outputSuffix);
        scorer.score(annotationsFile, outputScoresFile);
        return outputScoresFile;
    }
}