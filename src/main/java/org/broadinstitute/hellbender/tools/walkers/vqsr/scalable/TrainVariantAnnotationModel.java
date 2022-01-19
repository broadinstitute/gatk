package org.broadinstitute.hellbender.tools.walkers.vqsr.scalable;

import com.google.common.primitives.Doubles;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hdf5.HDF5File;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
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
public final class TrainVariantAnnotationModel extends CommandLineProgram {

    @Argument(
            fullName = "annotations-hdf5",
            doc = "HDF5 file containing annotations extracted with ExtractAnnotations.")
    private File annotationsHDF5File;

    @Argument(
            fullName = "python-script")
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

        IOUtils.canReadFile(pythonScriptFile);
        IOUtils.canReadFile(annotationsHDF5File);
        IOUtils.canReadFile(hyperparametersJSONFile);

        PythonScriptExecutor.checkPythonEnvironmentForPackage("sklearn");
        PythonScriptExecutor.checkPythonEnvironmentForPackage("dill");

        // TODO fail early for outputs
        final File outputTranchesFile = new File(outputPrefix + ".tranches.csv");
        final PrintStream tranchesStream;
        try {
            tranchesStream = new PrintStream(outputTranchesFile);
        } catch (final FileNotFoundException e) {
            throw new UserException.CouldNotCreateOutputFile(outputTranchesFile, e);
        }

        logger.info("Starting training...");

        final PythonScriptExecutor executor = new PythonScriptExecutor(true);
        final ProcessOutput pythonProcessOutput = executor.executeScriptAndGetOutput(
                pythonScriptFile.getAbsolutePath(),
                null,
                composePythonArguments(annotationsHDF5File, hyperparametersJSONFile, outputPrefix));

        if (pythonProcessOutput.getExitValue() != 0) {
            throw executor.getScriptException(executor.getExceptionMessageFromScriptError(pythonProcessOutput));
        }

        logger.info("Training complete.");

        final File outputScoresFile = new File(outputPrefix + ".scores.hdf5");
        try (final HDF5File outputScoresFileHDF5File = new HDF5File(outputScoresFile, HDF5File.OpenMode.READ_ONLY);
             final HDF5File inputAnnotationsHDF5File = new HDF5File(annotationsHDF5File, HDF5File.OpenMode.READ_ONLY)) {
            IOUtils.canReadFile(outputScoresFileHDF5File.getFile());
            final List<Double> scores = Doubles.asList(outputScoresFileHDF5File.readDoubleArray("/scores"));
            final List<Boolean> isTransition = Arrays.stream(inputAnnotationsHDF5File.readDoubleArray("/data/is_transition"))
                    .mapToObj(d -> (d == 1))
                    .collect(Collectors.toList());
            final List<Boolean> isTruth = Arrays.stream(inputAnnotationsHDF5File.readDoubleArray("/data/is_truth"))
                    .mapToObj(d -> (d == 1))
                    .collect(Collectors.toList());

            // Find the score cutoff values which correspond to the various tranches of calls requested by the user
            final int nCallsAtTruth = TruthSensitivityTranche.countCallsAtTruth(scores, isTruth, Double.NEGATIVE_INFINITY);
            final TruthSensitivityTranche.TruthSensitivityMetric metric = new TruthSensitivityTranche.TruthSensitivityMetric(nCallsAtTruth);
            final List<TruthSensitivityTranche> tranches = TruthSensitivityTranche.findTranches(scores, isTransition, isTruth, truthSensitivityTranches, metric, mode);
            tranchesStream.print(TruthSensitivityTranche.printHeader());
            tranchesStream.print(TruthSensitivityTranche.tranchesString(tranches));
        } catch (final RuntimeException exception) {
            throw new GATKException(String.format("Exception encountered during reading of scores from %s or annotations from %s: %s",
                    outputScoresFile.getAbsolutePath(), annotationsHDF5File.getAbsolutePath(), exception));
        }

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