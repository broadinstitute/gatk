package org.broadinstitute.hellbender.tools.walkers.vqsr.scalable;

import com.google.common.collect.Streams;
import htsjdk.variant.variantcontext.Allele;
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
import org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.modeling.VariantAnnotationsModel;
import org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.modeling.BGMMVariantAnnotationsScorer;
import org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.modeling.VariantAnnotationsScorer;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.python.PythonScriptExecutor;
import picard.cmdline.programgroups.VariantFilteringProgramGroup;

import java.io.File;
import java.io.IOException;
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
    
    private static final int DEFAULT_CHUNK_DIVISOR = 16;
    private static final int DEFAULT_MAXIMUM_CHUNK_SIZE = HDF5Utils.MAX_NUMBER_OF_VALUES_PER_HDF5_MATRIX / DEFAULT_CHUNK_DIVISOR;

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
            fullName = "ignore-variant-type")
    private boolean ignoreVariantType = false;

    /**
     * Add truth sensitivity slices through the call set at the given values. The default values are 100.0, 99.9, 99.0, and 90.0
     * which will result in 4 estimated tranches in the final call set: the full set of calls (100% sensitivity at the accessible
     * sites in the truth set), a 99.9% truth sensitivity tranche, along with progressively smaller tranches at 99% and 90%.
     * Note: You must pass in each tranche as a separate value (e.g. -tranche 100.0 -tranche 99.9).
     */
    @Argument(
            fullName = "truth-sensitivity-tranche",
            shortName = "tranche",
            doc = "The levels of truth sensitivity at which to slice the data. (in percent, that is 1.0 for 1 percent)",
            optional = true)
    private List<Double> truthSensitivityTranches = new ArrayList<>(Arrays.asList(100.0, 99.9, 99.0, 90.0));

    private ModelMode modelMode;

    @Override
    protected Object doWork() {

        IOUtils.canReadFile(inputAnnotationsFile);
        IOUtils.canReadFile(hyperparametersJSONFile);

        // TODO fail early for outputs, extract constants
        final File outputTrainingScoresFile = new File(outputPrefix + ".trainingScores.hdf5");
        final File outputTruthScoresFile = new File(outputPrefix + ".truthScores.hdf5");

        if (pythonScriptFile != null) {
            logger.info("Python script was provided, running in PYTHON mode...");
            modelMode = ModelMode.PYTHON;

            IOUtils.canReadFile(pythonScriptFile);
            PythonScriptExecutor.checkPythonEnvironmentForPackage("sklearn");
            PythonScriptExecutor.checkPythonEnvironmentForPackage("dill");
        } else {
            logger.info("Python script was not provided, running in BGMM mode...");
            modelMode = ModelMode.BGMM;
        }

        final VariantAnnotationsModel model;
        switch (modelMode) {
            case PYTHON:
//                model = new PythonVariantAnnotationsModel(hyperparametersJSONFile);
                model = new BGMMVariantAnnotationsModel(hyperparametersJSONFile);
                break;
            case BGMM:
                model = new BGMMVariantAnnotationsModel(hyperparametersJSONFile);
                break;
            default:
                throw new GATKException.ShouldNeverReachHereException("Unknown model mode.");
        }

        // handle SNP and INDEL here
        logger.info("Starting training...");

        // TODO could do subsetting out of memory by iterating over batches
        final List<String> annotationNames = LabeledVariantAnnotationsData.readAnnotationNames(inputAnnotationsFile);
        final double[][] allData = LabeledVariantAnnotationsData.readAnnotations(inputAnnotationsFile);
        final List<Boolean> isTraining = LabeledVariantAnnotationsData.readLabel(inputAnnotationsFile, LabeledVariantAnnotationsData.TRAINING_LABEL);
        final int numAllData = isTraining.size();
        final int numTraining = isTraining.stream().mapToInt(x -> x ? 1 : 0).sum();
        
        if (numTraining == 0) {
            throw new UserException.BadInput("No training sites were found in the provided annotations.");
        }
        
        // TODO clean up
        if (ignoreVariantType) {
            logger.info(String.format("Training type-agnostic model with %d training sites...", numTraining));
            final double[][] trainingData = IntStream.range(0, numAllData).boxed().filter(isTraining::get).map(i -> allData[i]).toArray(double[][]::new);
            final File trainingAnnotationsFile = IOUtils.createTempFile("training.annot", ".hdf5");

            try (final HDF5File trainingSNPAnnotationsHDF5File = new HDF5File(trainingAnnotationsFile, HDF5File.OpenMode.CREATE)) {
                trainingSNPAnnotationsHDF5File.makeStringArray("/annotations/names", annotationNames.toArray(new String[0]));
                HDF5Utils.writeChunkedDoubleMatrix(trainingSNPAnnotationsHDF5File, "/annotations", trainingData, DEFAULT_MAXIMUM_CHUNK_SIZE);
            } catch (final HDF5LibException exception) {
                throw new GATKException(String.format("Exception encountered during writing of annotations (%s). Output file at %s may be in a bad state.",
                        exception, trainingAnnotationsFile.getAbsolutePath()));
            }
            
            model.trainAndSerialize(inputAnnotationsFile, outputPrefix);
        } else {
            // SNP
            final List<Boolean> isSNP = LabeledVariantAnnotationsData.readLabel(inputAnnotationsFile, "snp");
            final List<Boolean> isTrainingSNP = Streams.zip(isTraining.stream(), isSNP.stream(), (a, b) -> a && b).collect(Collectors.toList());
            final int numTrainingSNPs = isTrainingSNP.stream().mapToInt(x -> x ? 1 : 0).sum();

            if (numTrainingSNPs > 0) {
                logger.info(String.format("Training SNP model with %d training SNPs...", numTrainingSNPs));
                final double[][] trainingSNPData = IntStream.range(0, numAllData).boxed().filter(isTrainingSNP::get).map(i -> allData[i]).toArray(double[][]::new);
                final File trainingSNPAnnotationsFile = IOUtils.createTempFile("training.snp.annot", ".hdf5");
                
                try (final HDF5File trainingSNPAnnotationsHDF5File = new HDF5File(trainingSNPAnnotationsFile, HDF5File.OpenMode.CREATE)) {
                    trainingSNPAnnotationsHDF5File.makeStringArray("/annotations/names", annotationNames.toArray(new String[0]));
                    HDF5Utils.writeChunkedDoubleMatrix(trainingSNPAnnotationsHDF5File, "/annotations", trainingSNPData, DEFAULT_MAXIMUM_CHUNK_SIZE);
                } catch (final HDF5LibException exception) {
                    throw new GATKException(String.format("Exception encountered during writing of annotations (%s). Output file at %s may be in a bad state.",
                            exception, trainingSNPAnnotationsFile.getAbsolutePath()));
                }

                model.trainAndSerialize(trainingSNPAnnotationsFile, outputPrefix + ".snp");
            }

            // INDEL
            final List<Boolean> isTrainingIndel = Streams.zip(isTraining.stream(), isSNP.stream(), (a, b) -> a && !b).collect(Collectors.toList());
            final int numTrainingIndels = isTrainingIndel.stream().mapToInt(x -> x ? 1 : 0).sum();
            
            if (numTrainingIndels > 0) {
                logger.info(String.format("Training Indel model with %d training indels...", numTrainingIndels));
                final double[][] trainingIndelData = IntStream.range(0, numAllData).boxed().filter(isTrainingIndel::get).map(i -> allData[i]).toArray(double[][]::new);
                final File trainingIndelAnnotationsFile = IOUtils.createTempFile("training.indel.annot", ".hdf5");

                try (final HDF5File trainingIndelAnnotationsHDF5File = new HDF5File(trainingIndelAnnotationsFile, HDF5File.OpenMode.CREATE)) {
                    trainingIndelAnnotationsHDF5File.makeStringArray("/annotations/names", annotationNames.toArray(new String[0]));
                    HDF5Utils.writeChunkedDoubleMatrix(trainingIndelAnnotationsHDF5File, "/annotations", trainingIndelData, DEFAULT_MAXIMUM_CHUNK_SIZE);
                } catch (final HDF5LibException exception) {
                    throw new GATKException(String.format("Exception encountered during writing of annotations (%s). Output file at %s may be in a bad state.",
                            exception, trainingIndelAnnotationsFile.getAbsolutePath()));
                }

                model.trainAndSerialize(trainingIndelAnnotationsFile, outputPrefix + ".indel");
            }
        }

//        final VariantAnnotationsScorer scorer;
//        switch (modelMode) {
//            case PYTHON:
////                scorer = new BGMMVariantAnnotationsScorer(hyperparametersJSONFile);
//                scorer = BGMMVariantAnnotationsScorer.deserialize(outputPrefix);
//                break;
//            case BGMM:
//                scorer = BGMMVariantAnnotationsScorer.deserialize(outputPrefix);
//                break;
//            default:
//                throw new GATKException.ShouldNeverReachHereException("Unknown model mode.");
//        }

//        // generate scores and write to HDF5
//        final double[] trainingScores = scorer.scoreSamples(preprocessedTrainingData);
//        VariantAnnotationUtils.writeScores(outputTrainingScoresFile, trainingScores);
//
//        final double[] truthScores = bgmm.scoreSamples(preprocessedTruthData);
//        VariantAnnotationUtils.writeScores(outputTruthScoresFile, truthScores);
//
//        if (pythonScriptFile != null) {
//
//
//            final PythonScriptExecutor executor = new PythonScriptExecutor(true);
//            final ProcessOutput pythonProcessOutput = executor.executeScriptAndGetOutput(
//                    pythonScriptFile.getAbsolutePath(),
//                    null,
//                    composePythonArguments(inputAnnotationsFile, hyperparametersJSONFile, outputPrefix));
//
//            if (pythonProcessOutput.getExitValue() != 0) {
//                throw executor.getScriptException(executor.getExceptionMessageFromScriptError(pythonProcessOutput));
//            }
//        } else {
//
//
//
//        }
//
//        VariantAnnotationUtils.writeTruthSensitivityTranches(
//                outputTranchesFile, outputTruthScoresFile, inputAnnotationsFile, truthSensitivityTranches, mode);

        logger.info(String.format("%s complete.", getClass().getSimpleName()));

        return null;
    }

    private static List<String> composePythonArguments(final File rawAnnotationsFile,
                                                       final File hyperparametersJSONFile,
                                                       final String outputPrefix) {
        try {
            return new ArrayList<>(Arrays.asList(
                    "--raw_annotations_file=" + rawAnnotationsFile.getCanonicalPath(),
                    "--hyperparameters_json_file=" + hyperparametersJSONFile.getCanonicalPath(),
                    "--output_prefix=" + outputPrefix));
        } catch (final IOException e) {
            throw new UserException.BadInput(String.format("Encountered exception resolving canonical file paths: %s", e));
        }
    }
}