package org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.modeling;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.python.PythonScriptExecutor;
import org.broadinstitute.hellbender.utils.runtime.ProcessOutput;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Given an HDF5 file containing annotations for a training set (in the format specified by
 * {@link VariantAnnotationsModel#trainAndSerialize}), a Python script containing modeling code,
 * and a JSON file containing hyperparameters, the {@link #trainAndSerialize} method can be used to train a model.
 *
 * The modeling script is expected to generate the file {outputPrefix}.scorer.pkl. This file should contain
 * a pickled Python lambda function to be used for generating scores from annotations in a subsequent test set.
 * The lambda should have the signature:
 *
 *      lambda test_annotation_names_i, test_X_ni
 *
 * Here, test_annotation_names_i is a numpy array of strings containing the annotation names, and
 * test X_ni is a numpy matrix of float-valued annotations, with dimensions (number of data points) x (number of annotations).
 * The lambda should check the test annotation names against the training annotation names and
 * then return a numpy array of float-valued scores with length given by the number of data points.
 *
 * See org/broadinstitute/hellbender/tools/walkers/vqsr/scalable/isolation-forest.py for an example implementation.
 */
public final class PythonSklearnVariantAnnotationsModel implements VariantAnnotationsModel {

    private final File pythonScriptFile;
    private final File hyperparametersJSONFile;

    public PythonSklearnVariantAnnotationsModel(final File pythonScriptFile,
                                                final File hyperparametersJSONFile) {
        this.pythonScriptFile = pythonScriptFile;
        this.hyperparametersJSONFile = hyperparametersJSONFile;
    }

    @Override
    public void trainAndSerialize(final File trainingAnnotationsFile,
                                  final String outputPrefix) {
        final PythonScriptExecutor executor = new PythonScriptExecutor(true);
        final ProcessOutput pythonProcessOutput = executor.executeScriptAndGetOutput(
                pythonScriptFile.getAbsolutePath(),
                null,
                composePythonArguments(trainingAnnotationsFile, hyperparametersJSONFile, outputPrefix));

        if (pythonProcessOutput.getExitValue() != 0) {
            throw executor.getScriptException(executor.getExceptionMessageFromScriptError(pythonProcessOutput));
        }
    }

    private static List<String> composePythonArguments(final File annotationsFile,
                                                       final File hyperparametersJSONFile,
                                                       final String outputPrefix) {
        try {
            return new ArrayList<>(Arrays.asList(
                    "--annotations_file=" + annotationsFile.getCanonicalPath(),
                    "--hyperparameters_json_file=" + hyperparametersJSONFile.getCanonicalPath(),
                    "--output_prefix=" + outputPrefix));
        } catch (final IOException e) {
            throw new UserException.BadInput(String.format("Encountered exception resolving canonical file paths: %s", e));
        }
    }
}