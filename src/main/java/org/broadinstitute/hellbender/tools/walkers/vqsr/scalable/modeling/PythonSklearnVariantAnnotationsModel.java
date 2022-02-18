package org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.modeling;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.python.PythonScriptExecutor;
import org.broadinstitute.hellbender.utils.runtime.ProcessOutput;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

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