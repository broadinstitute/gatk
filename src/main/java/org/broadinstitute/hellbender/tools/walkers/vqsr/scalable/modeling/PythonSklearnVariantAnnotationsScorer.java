package org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.modeling;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.python.PythonScriptExecutor;
import org.broadinstitute.hellbender.utils.runtime.ProcessOutput;

import java.io.File;
import java.io.IOException;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public final class PythonSklearnVariantAnnotationsScorer implements VariantAnnotationsScorer, Serializable {
    private static final long serialVersionUID = 1L;

    private final File pythonScriptFile;
    private final File scorerPklFile;

    public PythonSklearnVariantAnnotationsScorer(final File pythonScriptFile,
                                                 final File scorerPklFile) {
        this.pythonScriptFile = pythonScriptFile;
        this.scorerPklFile = scorerPklFile;
    }

    @Override
    public void score(final File inputAnnotationsFile,
                      final File outputScoresFile) {
        final PythonScriptExecutor executor = new PythonScriptExecutor(true);
        final ProcessOutput pythonProcessOutput = executor.executeScriptAndGetOutput(
                pythonScriptFile.getAbsolutePath(),
                null,
                composePythonArguments(inputAnnotationsFile, scorerPklFile, outputScoresFile));

        if (pythonProcessOutput.getExitValue() != 0) {
            throw executor.getScriptException(executor.getExceptionMessageFromScriptError(pythonProcessOutput));
        }
    }

    private static List<String> composePythonArguments(final File annotationsFile,
                                                       final File scorerPklFile,
                                                       final File outputScoresFile) {
        try {
            return new ArrayList<>(Arrays.asList(
                    "--annotations_file=" + annotationsFile.getCanonicalPath(),
                    "--scorer_pkl_file=" + scorerPklFile.getCanonicalPath(),
                    "--output_scores_file=" + outputScoresFile.getCanonicalPath()));
        } catch (final IOException e) {
            throw new UserException.BadInput(String.format("Encountered exception resolving canonical file paths: %s", e));
        }
    }
}
