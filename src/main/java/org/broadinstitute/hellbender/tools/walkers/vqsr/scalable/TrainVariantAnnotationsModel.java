package org.broadinstitute.hellbender.tools.walkers.vqsr.scalable;

import com.google.common.collect.Streams;
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
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.python.PythonScriptExecutor;
import picard.cmdline.programgroups.VariantFilteringProgramGroup;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

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

    @Argument(
            fullName = "annotations-hdf5",
            doc = "HDF5 file containing annotations extracted with ExtractAnnotations.")
    private File inputAnnotationsFile;

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
        final List<String> annotationNames = LabeledVariantAnnotationsData.readAnnotationNames(inputAnnotationsFile);
        final double[][] allAnnotations = LabeledVariantAnnotationsData.readAnnotations(inputAnnotationsFile);
        
        final List<Boolean> isTraining = LabeledVariantAnnotationsData.readLabel(inputAnnotationsFile, LabeledVariantAnnotationsData.TRAINING_LABEL);
        final List<Boolean> isTruth = LabeledVariantAnnotationsData.readLabel(inputAnnotationsFile, LabeledVariantAnnotationsData.TRUTH_LABEL);

        // TODO extract tags
        final List<Boolean> isSNP = LabeledVariantAnnotationsData.readLabel(inputAnnotationsFile, "snp");
        if (variantTypes.contains(VariantType.SNP)) {
            doModelingAndScoringWork(annotationNames, allAnnotations, isTraining, isTruth, isSNP, "SNP", ".snp");
        }
        if (variantTypes.contains(VariantType.INDEL)) {
            final List<Boolean> isIndel = isSNP.stream().map(x -> !x).collect(Collectors.toList());
            doModelingAndScoringWork(annotationNames, allAnnotations, isTraining, isTruth, isIndel, "INDEL", ".indel");
        }

        logger.info(String.format("%s complete.", getClass().getSimpleName()));

        return null;
    }

    private void doModelingAndScoringWork(final List<String> annotationNames,
                                          final double[][] allAnnotations,
                                          final List<Boolean> isTraining,
                                          final List<Boolean> isTruth,
                                          final List<Boolean> isVariantType,
                                          final String logMessageTag,
                                          final String outputPrefixTag) {
        final List<Boolean> isTrainingAndVariantType = Streams.zip(isTraining.stream(), isVariantType.stream(), (a, b) -> a && b).collect(Collectors.toList());
        final int numTrainingAndVariantType = numPassingFilter(isTrainingAndVariantType);

        if (numTrainingAndVariantType > 0) {
            logger.info(String.format("Training %s model with %d training sites x %d annotations %s...",
                    logMessageTag, numTrainingAndVariantType, annotationNames.size(), annotationNames));
            final File trainingAnnotationsFile = LabeledVariantAnnotationsData.subsetAnnotationsToTemporaryFile(annotationNames, allAnnotations, isTrainingAndVariantType);
            trainAndSerializeModel(trainingAnnotationsFile, outputPrefixTag);
            logger.info(String.format("Model trained and serialized with output prefix \"%s\".", outputPrefix + outputPrefixTag));

            final File trainingScoresFile = score(trainingAnnotationsFile, outputPrefixTag, TRAINING_SCORES_HDF5_SUFFIX);
            logger.info(String.format("Training scores written to %s.", trainingScoresFile.getAbsolutePath()));

            final List<Boolean> isTruthAndVariantType = Streams.zip(isTruth.stream(), isVariantType.stream(), (a, b) -> a && b).collect(Collectors.toList());
            final int numTruthAndVariantType = numPassingFilter(isTruthAndVariantType);
            if (numTruthAndVariantType > 0) {
                final File truthAnnotationsFile = LabeledVariantAnnotationsData.subsetAnnotationsToTemporaryFile(annotationNames, allAnnotations, isTruthAndVariantType);
                final File truthScoresFile = score(truthAnnotationsFile, outputPrefixTag, TRUTH_SCORES_HDF5_SUFFIX);
                logger.info(String.format("Truth scores written to %s.", truthScoresFile.getAbsolutePath()));
            }
        } else {
            throw new UserException.BadInput(String.format("Attempted to train %s model, " +
                    "but no suitable training sites were found in the provided annotations.", logMessageTag));
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
        scorer.scoreSamples(annotationsFile, outputScoresFile);
        return outputScoresFile;
    }
}