package org.broadinstitute.hellbender.tools.walkers.vqsr.scalable;

import com.google.common.primitives.Doubles;
import org.apache.commons.math3.linear.MatrixUtils;
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
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
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
public final class TrainVariantAnnotationModel extends CommandLineProgram {

    private static final String SCORER_SER_SUFFIX = ".scorer.ser";
    private static final String BGMM_HDF5_SUFFIX = ".bgmm.hdf5";

    @Argument(
            fullName = "annotations-hdf5",
            doc = "HDF5 file containing annotations extracted with ExtractAnnotations.")
    private File inputAnnotationsHDF5File;

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

    @Override
    protected Object doWork() {

        IOUtils.canReadFile(inputAnnotationsHDF5File);
        IOUtils.canReadFile(hyperparametersJSONFile);

        if (pythonScriptFile != null) {
            IOUtils.canReadFile(pythonScriptFile);
            PythonScriptExecutor.checkPythonEnvironmentForPackage("sklearn");
            PythonScriptExecutor.checkPythonEnvironmentForPackage("dill");
        }

        // TODO fail early for outputs
        final File outputTranchesFile = new File(outputPrefix + ".tranches.csv");
        final PrintStream tranchesStream;
        try {
            tranchesStream = new PrintStream(outputTranchesFile);
        } catch (final FileNotFoundException e) {
            throw new UserException.CouldNotCreateOutputFile(outputTranchesFile, e);
        }
        final File outputScoresFile = new File(outputPrefix + ".scores.hdf5");

        logger.info("Starting training...");

        if (pythonScriptFile != null) {
            logger.info("Python script was provided, running in PYTHON mode...");

            final PythonScriptExecutor executor = new PythonScriptExecutor(true);
            final ProcessOutput pythonProcessOutput = executor.executeScriptAndGetOutput(
                    pythonScriptFile.getAbsolutePath(),
                    null,
                    composePythonArguments(inputAnnotationsHDF5File, hyperparametersJSONFile, outputPrefix));

            if (pythonProcessOutput.getExitValue() != 0) {
                throw executor.getScriptException(executor.getExceptionMessageFromScriptError(pythonProcessOutput));
            }
        } else {
            logger.info("Python script was not provided, running in BGMM mode...");

            // load annotation training data
            logger.info("Reading annotations...");
            final double[][] data;
            try (final HDF5File annotationsHDF5File = new HDF5File(inputAnnotationsHDF5File, HDF5File.OpenMode.READ_ONLY)) {
                IOUtils.canReadFile(annotationsHDF5File.getFile());
                final List<Boolean> isTraining = readBooleanList(annotationsHDF5File, "/data/is_training");
                final double[][] allData = HDF5Utils.readChunkedDoubleMatrix(annotationsHDF5File, "/data/annotations");
                data = IntStream.range(0, isTraining.size()).boxed()
                        .filter(isTraining::get)
                        .map(i -> allData[i])
                        .toArray(double[][]::new);
            } catch (final HDF5LibException exception) {
                throw new GATKException(String.format("Exception encountered during reading of annotations from %s: %s",
                        inputAnnotationsHDF5File.getAbsolutePath(), exception));
            }
            final int nSamples = data.length;
            final int nFeatures = data[0].length;
            logger.info(String.format("Training BayesianGaussianMixtureModeller with %d training sites x %d annotations...", nSamples, nFeatures));

            // preprocess
            final BayesianGaussianMixtureUtils.Preprocesser preprocesser = new BayesianGaussianMixtureUtils.Preprocesser();
            final double[][] preprocessedData = preprocesser.fitTransform(data);

            // initialize BGMM

            final BayesianGaussianMixtureUtils.Hyperparameters hyperparameters =
                    BayesianGaussianMixtureUtils.Hyperparameters.readHyperparameters(hyperparametersJSONFile);
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

            // serialize preprocesser, BGMM, and scores
            // TODO fix up output paths and validation
            BayesianGaussianMixtureUtils.serialize(
                    new BayesianGaussianMixtureUtils.Scorer(preprocesser, bgmm),
                    new File(outputPrefix + SCORER_SER_SUFFIX));
            final BayesianGaussianMixtureModelPosterior fit = bgmm.getBestFit();
            fit.write(new File(outputPrefix + BGMM_HDF5_SUFFIX), "/bgmm");

            System.out.println("weights: " + fit.getWeights());
            System.out.println("meanPrecision: " + fit.getMeanPrecision());
            System.out.println("means: " + fit.getMeans());
            System.out.println("precisionsCholesky: " + fit.getPrecisionsCholesky());
            System.out.println("covariances: " + fit.getCovariances());
            System.out.println("degreesOfFreedom: " + fit.getDegreesOfFreedom());

            // write scores to HDF5
            final double[] scores = bgmm.scoreSamples(preprocessedData);
            try (final HDF5File outputScoresHDF5File = new HDF5File(outputScoresFile, HDF5File.OpenMode.CREATE)) {
                outputScoresHDF5File.makeDoubleArray("/data/scores", scores);
            } catch (final HDF5LibException exception) {
                throw new GATKException(String.format("Exception encountered during writing of scores (%s). Output file at %s may be in a bad state.",
                        exception, outputScoresFile.getAbsolutePath()));
            }
            logger.info(String.format("Scores written to %s.", outputScoresFile.getAbsolutePath()));
        }

        logger.info("Training complete.");

        // TODO some duplication of code here and in scoring tool
        try (final HDF5File outputScoresHDF5File = new HDF5File(outputScoresFile, HDF5File.OpenMode.READ_ONLY);
             final HDF5File inputAnnotationsHDF5File = new HDF5File(this.inputAnnotationsHDF5File, HDF5File.OpenMode.READ_ONLY)) {
            IOUtils.canReadFile(outputScoresHDF5File.getFile());
            IOUtils.canReadFile(inputAnnotationsHDF5File.getFile());
            final List<Double> scores = Doubles.asList(outputScoresHDF5File.readDoubleArray("/data/scores"));
            final List<Boolean> isBiallelicSNP = readBooleanList(inputAnnotationsHDF5File, "/data/is_biallelic_snp");
            final List<Boolean> isTransition = readBooleanList(inputAnnotationsHDF5File, "/data/is_transition");
            final List<Boolean> isTruth = readBooleanList(inputAnnotationsHDF5File, "/data/is_truth");

            // Find the score cutoff values which correspond to the various tranches of calls requested by the user
            final int nCallsAtTruth = TruthSensitivityTranche.countCallsAtTruth(scores, isTruth, Double.NEGATIVE_INFINITY);
            final TruthSensitivityTranche.TruthSensitivityMetric metric = new TruthSensitivityTranche.TruthSensitivityMetric(nCallsAtTruth);
            final List<TruthSensitivityTranche> tranches = TruthSensitivityTranche.findTranches(scores, isBiallelicSNP, isTransition, isTruth, truthSensitivityTranches, metric, mode);
            tranchesStream.print(TruthSensitivityTranche.printHeader());
            tranchesStream.print(TruthSensitivityTranche.tranchesString(tranches));
        } catch (final HDF5LibException exception) {
            throw new GATKException(String.format("Exception encountered during reading of scores from %s or annotations from %s: %s",
                    outputScoresFile.getAbsolutePath(), inputAnnotationsHDF5File.getAbsolutePath(), exception));
        }

        logger.info(String.format("%s complete.", getClass().getSimpleName()));

        return null;
    }

    private static List<Boolean> readBooleanList(final HDF5File annotationsHDF5File,
                                                 final String path) {
        return Arrays.stream(annotationsHDF5File.readDoubleArray(path))
                .mapToObj(d -> (d == 1))
                .collect(Collectors.toList());
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