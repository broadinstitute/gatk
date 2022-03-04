package org.broadinstitute.hellbender.tools.walkers.vqsr.scalable;

import com.google.common.collect.Streams;
import org.apache.commons.math3.stat.descriptive.rank.Percentile;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.data.LabeledVariantAnnotationsData;
import org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.data.VariantType;
import org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.modeling.BGMMVariantAnnotationsModel;
import org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.modeling.BGMMVariantAnnotationsScorer;
import org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.modeling.PythonSklearnVariantAnnotationsModel;
import org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.modeling.PythonSklearnVariantAnnotationsScorer;
import org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.modeling.VariantAnnotationsModel;
import org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.modeling.VariantAnnotationsScorer;
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

    enum ModelMode {
        PYTHON, BGMM
    }

    private static final String TRAINING_SCORES_HDF5_SUFFIX = ".trainingScores.hdf5";
    private static final String TRUTH_SCORES_HDF5_SUFFIX = ".truthScores.hdf5";
    private static final String UNLABELED_SCORES_HDF5_SUFFIX = ".unlabeledScores.hdf5";

    @Argument(
            fullName = "annotations-hdf5",
            doc = "HDF5 file containing annotations extracted with ExtractVariantAnnotations.")
    private File inputAnnotationsFile;

    @Argument(
            fullName = "unlabeled-annotations-hdf5",
            optional = true,
            doc = "HDF5 file containing annotations extracted with ExtractVariantAnnotations.")
    private File inputUnlabeledAnnotationsFile;

    @Argument(
            fullName = "python-script",
            optional = true)
    private File pythonScriptFile;

    @Argument(
            fullName = "hyperparameters-json",
            doc = "JSON file containing hyperparameters.")
    private File hyperparametersJSONFile;

    @Argument(
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "Output prefix.")
    private String outputPrefix;

    @Argument(
            fullName = "truth-sensitivity-threshold",
            doc = "",
            optional = true,
            minValue = 0.,
            maxValue = 1.) // TODO
    private Double truthSensitivityThreshold;

    @Argument(
            fullName = "mode",
            doc = "Variant types for which to train models. Duplicate values will be ignored.",
            minElements = 1,
            optional = true)
    public List<VariantType> variantTypes = new ArrayList<>(Arrays.asList(VariantType.SNP, VariantType.INDEL));

    private ModelMode modelMode;

    @Override
    protected Object doWork() {

        IOUtils.canReadFile(inputAnnotationsFile);
        IOUtils.canReadFile(hyperparametersJSONFile);

        if (inputUnlabeledAnnotationsFile != null) {
            IOUtils.canReadFile(inputUnlabeledAnnotationsFile);
        }

        // TODO test output and fail early

        if (pythonScriptFile != null) {
            logger.info("Python script was provided, running in PYTHON mode...");
            modelMode = ModelMode.PYTHON;

            IOUtils.canReadFile(pythonScriptFile);
            PythonScriptExecutor.checkPythonEnvironmentForPackage("argparse");
            PythonScriptExecutor.checkPythonEnvironmentForPackage("h5py");
            PythonScriptExecutor.checkPythonEnvironmentForPackage("numpy");
            PythonScriptExecutor.checkPythonEnvironmentForPackage("sklearn");
            PythonScriptExecutor.checkPythonEnvironmentForPackage("dill");
        } else {
            logger.info("Python script was not provided, running in BGMM mode...");
            modelMode = ModelMode.BGMM;
        }

        logger.info("Starting training...");

        // TODO could subset without allocating memory for allData by iterating over batches in inputAnnotationsFile
        final List<String> labeledAnnotationNames = LabeledVariantAnnotationsData.readAnnotationNames(inputAnnotationsFile);
        final double[][] labeledAnnotations = LabeledVariantAnnotationsData.readAnnotations(inputAnnotationsFile);
        
        final List<Boolean> labeledIsTraining = LabeledVariantAnnotationsData.readLabel(inputAnnotationsFile, LabeledVariantAnnotationsData.TRAINING_LABEL);
        final List<Boolean> labeledIsTruth = LabeledVariantAnnotationsData.readLabel(inputAnnotationsFile, LabeledVariantAnnotationsData.TRUTH_LABEL);

        // TODO extract tags
        final List<Boolean> labeledIsSNP = LabeledVariantAnnotationsData.readLabel(inputAnnotationsFile, "snp");
        if (variantTypes.contains(VariantType.SNP)) {
            doModelingAndScoringWork(labeledAnnotationNames, labeledAnnotations, labeledIsTraining, labeledIsTruth, labeledIsSNP, "SNP", ".snp");
        }
        if (variantTypes.contains(VariantType.INDEL)) {
            final List<Boolean> labeledIsIndel = labeledIsSNP.stream().map(x -> !x).collect(Collectors.toList());
            doModelingAndScoringWork(labeledAnnotationNames, labeledAnnotations, labeledIsTraining, labeledIsTruth, labeledIsIndel, "INDEL", ".indel");
        }

        // TODO extract, validate (including against truth-sensitivity-threshold and existence of truth scores)
        if (inputUnlabeledAnnotationsFile != null && truthSensitivityThreshold != null) {
            final List<String> unlabeledAnnotationNames = LabeledVariantAnnotationsData.readAnnotationNames(inputUnlabeledAnnotationsFile); // TODO validate
            final double[][] unlabeledAnnotations = LabeledVariantAnnotationsData.readAnnotations(inputUnlabeledAnnotationsFile);
            final List<Boolean> unlabeledIsSNP = LabeledVariantAnnotationsData.readLabel(inputUnlabeledAnnotationsFile, "snp");
            if (variantTypes.contains(VariantType.SNP)) {
                doNegativeModelingAndScoringWork(unlabeledAnnotationNames, unlabeledAnnotations, unlabeledIsSNP,
                        labeledAnnotations, labeledIsTraining, labeledIsTruth, labeledIsSNP, "SNP", ".snp");
            }
            if (variantTypes.contains(VariantType.INDEL)) {
                final List<Boolean> labeledIsIndel = labeledIsSNP.stream().map(x -> !x).collect(Collectors.toList());
                final List<Boolean> unlabeledIsIndel = unlabeledIsSNP.stream().map(x -> !x).collect(Collectors.toList());
                doNegativeModelingAndScoringWork(unlabeledAnnotationNames, unlabeledAnnotations, unlabeledIsIndel,
                        labeledAnnotations, labeledIsTraining, labeledIsTruth, labeledIsIndel,"INDEL", ".indel");
            }
        }

        logger.info(String.format("%s complete.", getClass().getSimpleName()));

        return null;
    }

    private void doModelingAndScoringWork(final List<String> annotationNames,
                                          final double[][] labeledAnnotations,
                                          final List<Boolean> isLabeledTraining,
                                          final List<Boolean> isLabeledTruth,
                                          final List<Boolean> isLabeledVariantType,
                                          final String logMessageTag,
                                          final String outputPrefixTag) {
        final List<Boolean> isLabeledTrainingAndVariantType = Streams.zip(isLabeledTraining.stream(), isLabeledVariantType.stream(), (a, b) -> a && b).collect(Collectors.toList());
        final int numLabeledTrainingAndVariantType = numPassingFilter(isLabeledTrainingAndVariantType);

        if (numLabeledTrainingAndVariantType > 0) {
            logger.info(String.format("Training %s model with %d training sites x %d annotations %s...",
                    logMessageTag, numLabeledTrainingAndVariantType, annotationNames.size(), annotationNames));
            final File labeledTrainingAndVariantTypeAnnotationsFile = LabeledVariantAnnotationsData.subsetAnnotationsToTemporaryFile(annotationNames, labeledAnnotations, isLabeledTrainingAndVariantType);
            trainAndSerializeModel(labeledTrainingAndVariantTypeAnnotationsFile, outputPrefixTag);
            logger.info(String.format("%s model trained and serialized with output prefix \"%s\".", logMessageTag, outputPrefix + outputPrefixTag));

            // TODO if BGMM, preprocess annotations and write to HDF5

            logger.info(String.format("Scoring %d %s training sites...", numLabeledTrainingAndVariantType, logMessageTag));
            final File labeledTrainingAndVariantTypeScoresFile = score(labeledTrainingAndVariantTypeAnnotationsFile, outputPrefixTag, TRAINING_SCORES_HDF5_SUFFIX);
            logger.info(String.format("%s training scores written to %s.", logMessageTag, labeledTrainingAndVariantTypeScoresFile.getAbsolutePath()));

            final List<Boolean> isLabeledTruthAndVariantType = Streams.zip(isLabeledTruth.stream(), isLabeledVariantType.stream(), (a, b) -> a && b).collect(Collectors.toList());
            final int numLabeledTruthAndVariantType = numPassingFilter(isLabeledTruthAndVariantType);
            if (numLabeledTruthAndVariantType > 0) {
                logger.info(String.format("Scoring %d %s truth sites...", numLabeledTruthAndVariantType, logMessageTag));
                final File labeledTruthAndVariantTypeAnnotationsFile = LabeledVariantAnnotationsData.subsetAnnotationsToTemporaryFile(annotationNames, labeledAnnotations, isLabeledTruthAndVariantType);
                final File labeledTruthAndVariantTypeScoresFile = score(labeledTruthAndVariantTypeAnnotationsFile, outputPrefixTag, TRUTH_SCORES_HDF5_SUFFIX);
                logger.info(String.format("%s truth scores written to %s.", logMessageTag, labeledTruthAndVariantTypeScoresFile.getAbsolutePath()));
            }
        } else {
            throw new UserException.BadInput(String.format("Attempted to train %s model, " +
                    "but no suitable training sites were found in the provided annotations.", logMessageTag));
        }
    }

    private void doNegativeModelingAndScoringWork(final List<String> annotationNames,
                                                  final double[][] unlabeledAnnotations,
                                                  final List<Boolean> isUnlabeledVariantType,
                                                  final double[][] labeledAnnotations,
                                                  final List<Boolean> isLabeledTraining,
                                                  final List<Boolean> isLabeledTruth,
                                                  final List<Boolean> isLabeledVariantType,
                                                  final String logMessageTag,
                                                  final String outputPrefixTag) {
        final int numVariantType = numPassingFilter(isUnlabeledVariantType);

        if (numVariantType > 0) {
            final File labeledTruthAndVariantTypeScoresFile = new File(outputPrefix + outputPrefixTag + TRUTH_SCORES_HDF5_SUFFIX); // produced by doModelingAndScoringWork TODO output a copy of these?
            final double[] labeledTruthAndVariantTypeScores = VariantAnnotationsScorer.readScores(labeledTruthAndVariantTypeScoresFile);
            final double scoreThreshold = new Percentile(100. * (1. - truthSensitivityThreshold)).evaluate(labeledTruthAndVariantTypeScores);
            logger.info(String.format("Using %s score threshold of %.4f corresponding to specified truth-sensitivity threshold of %.4f ...",
                    logMessageTag, scoreThreshold, truthSensitivityThreshold));

            final File labeledTrainingAndVariantTypeScoresFile = new File(outputPrefix + outputPrefixTag + TRAINING_SCORES_HDF5_SUFFIX); // produced by doModelingAndScoringWork TODO output a copy of these?
            final double[] labeledTrainingAndVariantTypeScores = VariantAnnotationsScorer.readScores(labeledTrainingAndVariantTypeScoresFile);
            final List<Boolean> isNegativeTrainingFromLabeledTrainingAndVariantType = Arrays.stream(labeledTrainingAndVariantTypeScores).boxed().map(s -> s < scoreThreshold).collect(Collectors.toList());
            final int numNegativeTrainingFromLabeledTrainingAndVariantType = numPassingFilter(isNegativeTrainingFromLabeledTrainingAndVariantType);
            logger.info(String.format("Selected %d labeled %s sites below score threshold of %.4f for negative-model training...",
                    numNegativeTrainingFromLabeledTrainingAndVariantType, logMessageTag, scoreThreshold));

            logger.info(String.format("Scoring %d unlabeled %s sites...", numVariantType, logMessageTag));
            final File unlabeledVariantTypeAnnotationsFile = LabeledVariantAnnotationsData.subsetAnnotationsToTemporaryFile(annotationNames, unlabeledAnnotations, isUnlabeledVariantType);
            final File unlabeledVariantTypeScoresFile = score(unlabeledVariantTypeAnnotationsFile, outputPrefixTag, UNLABELED_SCORES_HDF5_SUFFIX);
            final double[] unlabeledVariantTypeScores = VariantAnnotationsScorer.readScores(unlabeledVariantTypeScoresFile);
            final List<Boolean> isNegativeTrainingFromUnlabeledVariantType = Arrays.stream(unlabeledVariantTypeScores).boxed().map(s -> s < scoreThreshold).collect(Collectors.toList()); // length matches unlabeledAnnotationsFile
            final int numNegativeTrainingFromUnlabeledVariantType = numPassingFilter(isNegativeTrainingFromUnlabeledVariantType);
            logger.info(String.format("Selected %d unlabeled %s sites below score threshold of %.4f for negative-model training...",
                    numNegativeTrainingFromUnlabeledVariantType, logMessageTag, scoreThreshold));

            // combine labeled and unlabeled negative training data
            final File negativeTrainingFromLabeledTrainingAndVariantTypeAnnotationsFile =
                    LabeledVariantAnnotationsData.subsetAnnotationsToTemporaryFile(annotationNames, labeledAnnotations, isNegativeTrainingFromLabeledTrainingAndVariantType);
            final double[][] negativeTrainingFromLabeledTrainingAndVariantTypeAnnotations = LabeledVariantAnnotationsData.readAnnotations(negativeTrainingFromLabeledTrainingAndVariantTypeAnnotationsFile);

            final File negativeTrainingFromUnlabeledVariantTypeAnnotationsFile =
                    LabeledVariantAnnotationsData.subsetAnnotationsToTemporaryFile(annotationNames, unlabeledAnnotations, isNegativeTrainingFromUnlabeledVariantType);
            final double[][] negativeTrainingFromUnlabeledVariantTypeAnnotations = LabeledVariantAnnotationsData.readAnnotations(negativeTrainingFromUnlabeledVariantTypeAnnotationsFile);

            final double[][] negativeTrainingAnnotations = Streams.concat(
                    Arrays.stream(negativeTrainingFromLabeledTrainingAndVariantTypeAnnotations),
                    Arrays.stream(negativeTrainingFromUnlabeledVariantTypeAnnotations)).toArray(double[][]::new);
            final int numNegativeTraining = negativeTrainingAnnotations.length;
            final File negativeTrainingAnnotationsFile = LabeledVariantAnnotationsData.subsetAnnotationsToTemporaryFile(
                    annotationNames, negativeTrainingAnnotations, Collections.nCopies(numNegativeTraining, true));

            logger.info(String.format("Training %s negative model with %d negative-training sites x %d annotations %s...",
                    logMessageTag, numNegativeTraining, annotationNames.size(), annotationNames));
            trainAndSerializeModel(negativeTrainingAnnotationsFile, outputPrefixTag + ".negative");
            logger.info(String.format("%s negative model trained and serialized with output prefix \"%s\".", logMessageTag, outputPrefix + outputPrefixTag + ".negative"));

//            final List<Boolean> isLabeledTrainingAndVariantType = Streams.zip(isLabeledTraining.stream(), isLabeledVariantType.stream(), (a, b) -> a && b).collect(Collectors.toList());
//            final int numLabeledTrainingAndVariantType = numPassingFilter(isLabeledTrainingAndVariantType);
//            final File labeledTrainingAnnotationsFile = LabeledVariantAnnotationsData.subsetAnnotationsToTemporaryFile(annotationNames, labeledAnnotations, isLabeledTrainingAndVariantType);
//            logger.info(String.format("Re-scoring %d %s training sites...", numLabeledTrainingAndVariantType, logMessageTag));
//            final File labeledTrainingScoresFile = positiveNegativeScore(labeledTrainingAnnotationsFile, outputPrefixTag, TRAINING_SCORES_HDF5_SUFFIX);
//            logger.info(String.format("Training scores written to %s.", labeledTrainingScoresFile.getAbsolutePath()));

            final List<Boolean> isLabeledTruthAndVariantType = Streams.zip(isLabeledTruth.stream(), isLabeledVariantType.stream(), (a, b) -> a && b).collect(Collectors.toList());
            final int numLabeledTruthAndVariantType = numPassingFilter(isLabeledTruthAndVariantType);
            if (numLabeledTruthAndVariantType > 0) {
                logger.info(String.format("Re-scoring %d %s truth sites...", numLabeledTruthAndVariantType, logMessageTag));
                final File labeledTruthAnnotationsFile = LabeledVariantAnnotationsData.subsetAnnotationsToTemporaryFile(annotationNames, labeledAnnotations, isLabeledTruthAndVariantType);
                final File labeledTruthScoresFile = positiveNegativeScore(labeledTruthAnnotationsFile, outputPrefixTag, TRUTH_SCORES_HDF5_SUFFIX);
                logger.info(String.format("Truth scores written to %s.", labeledTruthScoresFile.getAbsolutePath()));
            }
        } else {
            throw new UserException.BadInput(String.format("Attempted to train %s negative model, " +
                    "but no suitable sites were found in the provided annotations.", logMessageTag));
        }
    }

    private static int numPassingFilter(List<Boolean> isPassing) {
        return isPassing.stream().mapToInt(x -> x ? 1 : 0).sum();
    }

    private void trainAndSerializeModel(final File trainingAnnotationsFile,
                                        final String outputPrefixTag) {
        final VariantAnnotationsModel model;
        switch (modelMode) {
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
        switch (modelMode) {
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
        switch (modelMode) {
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