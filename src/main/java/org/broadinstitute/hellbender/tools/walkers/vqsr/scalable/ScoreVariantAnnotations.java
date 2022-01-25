package org.broadinstitute.hellbender.tools.walkers.vqsr.scalable;

import com.google.common.primitives.Doubles;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hdf5.HDF5File;
import org.broadinstitute.hdf5.HDF5LibException;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
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
public class ScoreVariantAnnotations extends VariantAnnotationWalker {

    private static final String SCORER_PKL_SUFFIX = ".scorer.pkl";
    private static final String SCORER_SER_SUFFIX = ".scorer.ser";
    private static final String SCORES_HDF5_SUFFIX = ".scores.hdf5";
    private static final String RECALIBRATION_VCF_SUFFIX = ".recal.vcf";

    @Argument(
            fullName = "python-script",
            optional = true)
    private File pythonScriptFile;

    @Argument(
            fullName = "model-prefix")
    private String modelPrefix;

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

    private File inputScorerPklFile;
    private File outputScoresHDF5File;
    private PrintStream tranchesStream;

    @Override
    public String getVCFSuffix() {
        return RECALIBRATION_VCF_SUFFIX;
    }

    @Override
    public void beforeOnTraversalStart() {
        isExtractTrainingAndTruthOnly = false;

        if (pythonScriptFile != null) {
            logger.info("Python script was provided, running in PYTHON mode...");

            inputScorerPklFile = new File(modelPrefix + SCORER_PKL_SUFFIX);

            IOUtils.canReadFile(pythonScriptFile);
            IOUtils.canReadFile(inputScorerPklFile);

            PythonScriptExecutor.checkPythonEnvironmentForPackage("sklearn");
            PythonScriptExecutor.checkPythonEnvironmentForPackage("dill");
        } else {
            // TODO validate BGMM model inputs
        }

        outputScoresHDF5File = new File(outputPrefix + SCORES_HDF5_SUFFIX);

        for (final File outputFile : Arrays.asList(outputScoresHDF5File)) {
            if ((outputFile.exists() && !outputFile.canWrite()) ||
                    (!outputFile.exists() && !outputFile.getAbsoluteFile().getParentFile().canWrite())) {
                throw new UserException(String.format("Cannot create output file at %s.", outputFile));
            }
        }

        final File outputTranchesFile = new File(outputPrefix + ".tranches.csv");
        try {
            tranchesStream = new PrintStream(outputTranchesFile);
        } catch (final FileNotFoundException e) {
            throw new UserException.CouldNotCreateOutputFile(outputTranchesFile, e);
        }
    }

    @Override
    public void afterTraversalSuccess() {

        logger.info(String.format("Extracted annotations for %s total variants.", dataManager.getData().size()));

        logger.info("Writing annotations...");
        writeAnnotationsHDF5();

        logger.info("Scoring...");
        final double[] scores;
        if (pythonScriptFile != null) {
            final PythonScriptExecutor executor = new PythonScriptExecutor(true);
            final ProcessOutput pythonProcessOutput = executor.executeScriptAndGetOutput(
                    pythonScriptFile.getAbsolutePath(),
                    null,
                    composePythonArguments(outputAnnotationsHDF5File, inputScorerPklFile, outputScoresHDF5File));

            if (pythonProcessOutput.getExitValue() != 0) {
                throw executor.getScriptException(executor.getExceptionMessageFromScriptError(pythonProcessOutput));
            }

            try (final HDF5File outputScoresFileHDF5File = new HDF5File(outputScoresHDF5File, HDF5File.OpenMode.READ_ONLY)) {
                IOUtils.canReadFile(outputScoresFileHDF5File.getFile());
                scores = outputScoresFileHDF5File.readDoubleArray("/data/scores");
                VariantDataManager.setScores(dataManager.getData(), scores);
            } catch (final HDF5LibException exception) {
                throw new GATKException(String.format("Exception encountered during reading of scores from %s: %s",
                        outputScoresHDF5File.getAbsolutePath(), exception));
            }
        } else {
            final BayesianGaussianMixtureUtils.Scorer scorer = BayesianGaussianMixtureUtils.deserialize(
                    new File(modelPrefix + SCORER_SER_SUFFIX), // TODO clean up
                    BayesianGaussianMixtureUtils.Scorer.class);
            final double[][] data = dataManager.getData().stream().map(vd -> vd.annotations).toArray(double[][]::new);
            scores = scorer.preprocessAndScoreSamples(data);
        }
        logger.info("Scoring complete.");

        logger.info("Writing VCF...");
        writeVCF(false, false,true);

        // TODO we could just get all this stuff from the VariantDataManager, but let's clean that up later
        // TODO some duplication of code here and in training tool
        try (final HDF5File annotationsHDF5File = new HDF5File(outputAnnotationsHDF5File, HDF5File.OpenMode.READ_ONLY)) {
            IOUtils.canReadFile(annotationsHDF5File.getFile());
            final List<Boolean> isBiallelicSNP = readBooleanList(annotationsHDF5File, "/data/is_biallelic_snp");
            final List<Boolean> isTransition = readBooleanList(annotationsHDF5File, "/data/is_transition");
            final List<Boolean> isTruth = readBooleanList(annotationsHDF5File, "/data/is_truth");

            // Find the score cutoff values which correspond to the various tranches of calls requested by the user
            final int nCallsAtTruth = TruthSensitivityTranche.countCallsAtTruth(Doubles.asList(scores), isTruth, Double.NEGATIVE_INFINITY);
            final TruthSensitivityTranche.TruthSensitivityMetric metric = new TruthSensitivityTranche.TruthSensitivityMetric(nCallsAtTruth);
            final List<TruthSensitivityTranche> tranches = TruthSensitivityTranche.findTranches(Doubles.asList(scores), isBiallelicSNP, isTransition, isTruth, truthSensitivityTranches, metric, mode);
            tranchesStream.print(TruthSensitivityTranche.printHeader());
            tranchesStream.print(TruthSensitivityTranche.tranchesString(tranches));
        } catch (final HDF5LibException exception) {
            throw new GATKException(String.format("Exception encountered during reading of annotations from %s: %s",
                    outputAnnotationsHDF5File.getAbsolutePath(), exception));
        }

        logger.info(String.format("%s complete.", getClass().getSimpleName()));
    }

    private static List<Boolean> readBooleanList(final HDF5File annotationsHDF5File,
                                                 final String path) {
        return Arrays.stream(annotationsHDF5File.readDoubleArray(path))
                .mapToObj(d -> (d == 1))
                .collect(Collectors.toList());
    }

    private static List<String> composePythonArguments(final File rawAnnotationsFile,
                                                       final File scorerPklFile,
                                                       final File outputScoresFile) {
        try {
            return new ArrayList<>(Arrays.asList(
                    "--raw_annotations_file=" + rawAnnotationsFile.getCanonicalPath(),
                    "--scorer_pkl_file=" + scorerPklFile.getCanonicalPath(),
                    "--output_scores_file=" + outputScoresFile.getCanonicalPath()));
        } catch (final IOException e) {
            throw new UserException.BadInput(String.format("Encountered exception resolving canonical file paths: %s", e));
        }
    }
}