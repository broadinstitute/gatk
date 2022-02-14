package org.broadinstitute.hellbender.tools.walkers.vqsr.scalable;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.data.VariantAnnotationUtils;
import org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.data.VariantType;
import org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.model.BGMMVariantAnnotationsModeller;
import org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.model.VariantAnnotationsModeller;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.python.PythonScriptExecutor;
import org.broadinstitute.hellbender.utils.runtime.ProcessOutput;
import picard.cmdline.programgroups.VariantFilteringProgramGroup;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

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

    static enum ModelMode {
        PYTHON, BGMM
    }

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
    public VariantType mode = VariantType.SNP;

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

        final VariantAnnotationsModeller model;
        switch (modelMode) {
            case PYTHON:
                model = new PythonVariantAnnotationsModel(hyperparametersJSONFile);
                break;
            case BGMM:
                model = new BGMMVariantAnnotationsModeller(hyperparametersJSONFile, inputAnnotationsFile);
                break;
            default:
                throw new GATKException.ShouldNeverReachHereException("Unknown model mode.");
        }

        // handle SNP and INDEL here
        logger.info("Starting training...");
        model.fitAndWriteTrainingAndTruthScoresAndSerialize(inputAnnotationsFile, outputPrefix);

        // log output paths

        // serialize calibrated scorer if truth scores available

        // generate scores and write to HDF5
        final double[] trainingScores = bgmm.scoreSamples(preprocessedTrainingData);
        VariantAnnotationUtils.writeScores(outputTrainingScoresFile, trainingScores);

        final double[] truthScores = bgmm.scoreSamples(preprocessedTruthData);
        VariantAnnotationUtils.writeScores(outputTruthScoresFile, truthScores);

        if (pythonScriptFile != null) {


            final PythonScriptExecutor executor = new PythonScriptExecutor(true);
            final ProcessOutput pythonProcessOutput = executor.executeScriptAndGetOutput(
                    pythonScriptFile.getAbsolutePath(),
                    null,
                    composePythonArguments(inputAnnotationsFile, hyperparametersJSONFile, outputPrefix));

            if (pythonProcessOutput.getExitValue() != 0) {
                throw executor.getScriptException(executor.getExceptionMessageFromScriptError(pythonProcessOutput));
            }
        } else {



        }

        VariantAnnotationUtils.writeTruthSensitivityTranches(
                outputTranchesFile, outputTruthScoresFile, inputAnnotationsFile, truthSensitivityTranches, mode);

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