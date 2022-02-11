package org.broadinstitute.hellbender.tools.walkers.vqsr.scalable;

import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hdf5.HDF5File;
import org.broadinstitute.hdf5.HDF5LibException;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.utils.HDF5Utils;
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
    private File outputScoresFile;
    private File outputTranchesFile;

    @Override
    public String getVCFSuffix() {
        return RECALIBRATION_VCF_SUFFIX;
    }

    @Override
    public void beforeOnTraversalStart() {
        isExtractAll = true;

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

        outputScoresFile = new File(outputPrefix + SCORES_HDF5_SUFFIX);
        outputTranchesFile = new File(outputPrefix + ".tranches.csv");

        for (final File outputFile : Arrays.asList(outputScoresFile, outputTranchesFile)) {
            if ((outputFile.exists() && !outputFile.canWrite()) ||
                    (!outputFile.exists() && !outputFile.getAbsoluteFile().getParentFile().canWrite())) {
                throw new UserException(String.format("Cannot create output file at %s.", outputFile));
            }
        }
    }

    @Override
    public void afterTraversalSuccess() {

        logger.info(String.format("Extracted annotations for %s total variants.", data.getData().size()));

        logger.info("Writing annotations...");
        writeAnnotationsHDF5();

        logger.info("Scoring...");
        final double[] scores;
        if (pythonScriptFile != null) {
            final PythonScriptExecutor executor = new PythonScriptExecutor(true);
            final ProcessOutput pythonProcessOutput = executor.executeScriptAndGetOutput(
                    pythonScriptFile.getAbsolutePath(),
                    null,
                    composePythonArguments(outputAnnotationsFile, inputScorerPklFile, outputScoresFile));

            if (pythonProcessOutput.getExitValue() != 0) {
                throw executor.getScriptException(executor.getExceptionMessageFromScriptError(pythonProcessOutput));
            }
            scores = VariantAnnotationUtils.readScores(outputScoresFile);
        } else {
            final VariantAnnotationUtils.Scorer scorer = VariantAnnotationUtils.deserialize(
                    new File(modelPrefix + SCORER_SER_SUFFIX), // TODO clean up
                    VariantAnnotationUtils.Scorer.class);
            final double[][] data = this.data.getData().stream().map(vd -> vd.annotations).toArray(double[][]::new);
            final Pair<double[][], double[]> preprocessedDataAndScores = scorer.preprocessAndScoreSamples(data);
            final double[][] preprocessedData = preprocessedDataAndScores.getLeft();
            scores = preprocessedDataAndScores.getRight();

            // write preprocessed annotations
            // TODO clean this up
            final List<String> annotationNames = this.data.getSortedAnnotationKeys();

            final File outputPreprocessedAnnotationsFile = new File(outputPrefix + ".annot.pre.hdf5");
            try (final HDF5File hdf5File = new HDF5File(outputPreprocessedAnnotationsFile, HDF5File.OpenMode.CREATE)) { // TODO allow appending
                IOUtils.canReadFile(hdf5File.getFile());

                hdf5File.makeStringArray("/data/annotation_names", annotationNames.toArray(new String[0]));
                HDF5Utils.writeChunkedDoubleMatrix(hdf5File, "/data/annotations", preprocessedData, maximumChunkSize);
                hdf5File.makeDoubleArray("/data/is_training", this.data.getData().stream().mapToDouble(x -> x.labels.contains("training") ? 1 : 0).toArray());
            } catch (final HDF5LibException exception) {
                throw new GATKException(String.format("Exception encountered during writing of preprocessed annotations (%s). Output file at %s may be in a bad state.",
                        exception, outputPreprocessedAnnotationsFile.getAbsolutePath()));
            }
            logger.info(String.format("Preprocessed annotations written to %s.", outputPreprocessedAnnotationsFile.getAbsolutePath()));
        }
        VariantDataCollection.setScores(data.getData(), scores);
        logger.info("Scoring complete.");

        logger.info("Writing VCF...");
        writeVCF(false, false,true);

        VariantAnnotationUtils.writeTruthSensitivityTranches(
                outputTranchesFile, outputScoresFile, outputAnnotationsFile, truthSensitivityTranches, mode);

        logger.info(String.format("%s complete.", getClass().getSimpleName()));
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