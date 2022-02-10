package org.broadinstitute.hellbender.tools.walkers.vqsr.scalable;

import org.apache.commons.math3.linear.MatrixUtils;
import org.broadinstitute.barclay.argparser.Advanced;
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
import org.broadinstitute.hellbender.utils.clustering.BayesianGaussianMixtureModelPosterior;
import org.broadinstitute.hellbender.utils.clustering.BayesianGaussianMixtureModeller;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.python.PythonScriptExecutor;
import org.broadinstitute.hellbender.utils.runtime.ProcessOutput;
import picard.cmdline.programgroups.VariantFilteringProgramGroup;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
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
public final class TrainVariantAnnotationModel extends CommandLineProgram {

    private static final String SCORER_SER_SUFFIX = ".scorer.ser";
    private static final String BGMM_HDF5_SUFFIX = ".bgmm.hdf5";

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

    // TODO get rid of this argument, it's only needed for tranches
    @Argument(
            fullName = "mode",
            shortName = "mode",
            doc = "Variant type to extract")
    public VariantTypeMode mode = VariantTypeMode.SNP;

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

    @Advanced
    @Argument(
            doc = "Maximum HDF5 matrix chunk size.  Large matrices written to HDF5 are chunked into equally sized " +
                    "subsets of rows (plus a subset containing the remainder, if necessary) to avoid a hard limit in " +
                    "Java HDF5 on the number of elements in a matrix.  However, since a single row is not allowed to " +
                    "be split across multiple chunks, the number of columns must be less than the maximum number of " +
                    "values in each chunk.  Decreasing this number will reduce heap usage when writing chunks.",
            fullName = "maximum-chunk-size",
            minValue = 1,
            maxValue = HDF5Utils.MAX_NUMBER_OF_VALUES_PER_HDF5_MATRIX,
            optional = true
    )
    private int maximumChunkSize = VariantAnnotationWalker.DEFAULT_MAXIMUM_CHUNK_SIZE;

    @Override
    protected Object doWork() {

        IOUtils.canReadFile(inputAnnotationsFile);
        IOUtils.canReadFile(hyperparametersJSONFile);

        if (pythonScriptFile != null) {
            IOUtils.canReadFile(pythonScriptFile);
            PythonScriptExecutor.checkPythonEnvironmentForPackage("sklearn");
            PythonScriptExecutor.checkPythonEnvironmentForPackage("dill");
        }

        // TODO fail early for outputs
        final File outputTranchesFile = new File(outputPrefix + ".tranches.csv");
        final File outputScoresFile = new File(outputPrefix + ".scores.hdf5");

        logger.info("Starting training...");

        if (pythonScriptFile != null) {
            logger.info("Python script was provided, running in PYTHON mode...");

            final PythonScriptExecutor executor = new PythonScriptExecutor(true);
            final ProcessOutput pythonProcessOutput = executor.executeScriptAndGetOutput(
                    pythonScriptFile.getAbsolutePath(),
                    null,
                    composePythonArguments(inputAnnotationsFile, hyperparametersJSONFile, outputPrefix));

            if (pythonProcessOutput.getExitValue() != 0) {
                throw executor.getScriptException(executor.getExceptionMessageFromScriptError(pythonProcessOutput));
            }
        } else {
            logger.info("Python script was not provided, running in BGMM mode...");

            // load annotation training data
            logger.info("Reading annotations...");
            final double[][] allData = VariantDataCollection.readData(inputAnnotationsFile);
            final String[] annotationNames = VariantDataCollection.readAnnotationNames(inputAnnotationsFile);
            final List<Boolean> isTraining = VariantDataCollection.readLabel(inputAnnotationsFile, VariantDataCollection.TRAINING_LABEL);
            final double[][] data = IntStream.range(0, isTraining.size()).boxed()
                    .filter(isTraining::get)
                    .map(i -> allData[i])
                    .toArray(double[][]::new);
            final int nSamples = data.length;
            final int nFeatures = data[0].length;
            logger.info(String.format("Training BayesianGaussianMixtureModeller with %d training sites x %d annotations...", nSamples, nFeatures));

            // preprocess
            final VariantAnnotationUtils.Preprocesser preprocesser = new VariantAnnotationUtils.Preprocesser();
            final double[][] preprocessedData = preprocesser.fitTransform(data);

            // write preprocessed annotations
            // TODO clean this up
            final File outputPreprocessedAnnotationsFile = new File(outputPrefix + ".annot.pre.hdf5");
            try (final HDF5File hdf5File = new HDF5File(outputPreprocessedAnnotationsFile, HDF5File.OpenMode.CREATE)) { // TODO allow appending
                IOUtils.canReadFile(hdf5File.getFile());

                hdf5File.makeStringArray("/data/annotation_names", annotationNames);
                HDF5Utils.writeChunkedDoubleMatrix(hdf5File, "/data/annotations", preprocessedData, maximumChunkSize);
                hdf5File.makeDoubleArray("/data/is_training", isTraining.stream().mapToDouble(x -> x ? 1 : 0).toArray());
            } catch (final HDF5LibException exception) {
                throw new GATKException(String.format("Exception encountered during writing of preprocessed annotations (%s). Output file at %s may be in a bad state.",
                        exception, outputPreprocessedAnnotationsFile.getAbsolutePath()));
            }
            logger.info(String.format("Preprocessed annotations written to %s.", outputPreprocessedAnnotationsFile.getAbsolutePath()));

            // BGMM
            final VariantAnnotationUtils.Hyperparameters hyperparameters =
                    VariantAnnotationUtils.Hyperparameters.readHyperparameters(hyperparametersJSONFile);
            final double[][] covariancePrior = MatrixUtils.createRealIdentityMatrix(nFeatures).getData();
            final BayesianGaussianMixtureModeller bgmm = new BayesianGaussianMixtureModeller.Builder()
                    .nComponents(hyperparameters.nComponents)
                    .tol(hyperparameters.tol)
                    .regCovar(hyperparameters.regCovar)
                    .maxIter(hyperparameters.maxIter)
                    .nInit(hyperparameters.nInit)
                    .initMethod(hyperparameters.getInitMethod())
                    .weightConcentrationPrior(hyperparameters.weightConcentrationPrior)
                    .meanPrecisionPrior(hyperparameters.meanPrecisionPrior)
                    .meanPrior(
                            hyperparameters.meanPrior != null
                                    ? Collections.nCopies(nFeatures, hyperparameters.meanPrior).stream().mapToDouble(x -> x).toArray()
                                    : null)
                    .degreesOfFreedomPrior(hyperparameters.degreesOfFreedomPrior)
                    .covariancePrior(covariancePrior)
                    .seed(hyperparameters.seed)
                    .warmStart(hyperparameters.warmStart)
                    .verboseInterval(hyperparameters.verboseInterval)
                    .build();
            bgmm.fit(preprocessedData);

            // serialize scorer = preprocesser + BGMM
            // TODO fix up output paths and validation
            VariantAnnotationUtils.serialize(
                    new File(outputPrefix + SCORER_SER_SUFFIX), new VariantAnnotationUtils.Scorer(preprocesser, bgmm)
            );
            final BayesianGaussianMixtureModelPosterior fit = bgmm.getBestFit();
            fit.write(new File(outputPrefix + BGMM_HDF5_SUFFIX), "/bgmm");

            // generate scores and write to HDF5
            final double[] scores = bgmm.scoreSamples(preprocessedData);
            VariantAnnotationUtils.writeScores(outputScoresFile, scores);
        }

        logger.info("Training complete.");

        VariantAnnotationUtils.writeTruthSensitivityTranches(
                outputTranchesFile, outputScoresFile, inputAnnotationsFile, truthSensitivityTranches, mode);

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