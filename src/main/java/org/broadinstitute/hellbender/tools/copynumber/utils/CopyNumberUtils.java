package org.broadinstitute.hellbender.tools.copynumber.utils;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.python.PythonScriptExecutor;

import java.io.File;
import java.io.IOException;

public final class CopyNumberUtils {
    private CopyNumberUtils() {}

    /**
     * File paths that are passed to {@link PythonScriptExecutor} must be canonical (rather than absolute).
     * See https://github.com/broadinstitute/gatk/issues/4724.
     */
    public static String getCanonicalPath(final File file) {
        Utils.nonNull(file);
        try {
            return file.getCanonicalPath();
        } catch (final IOException e) {
            throw new UserException.BadInput(String.format("Could not resolve a canonical file path: %s", file));
        }
    }

    /**
     * File paths that are passed to {@link PythonScriptExecutor} must be canonical (rather than absolute).
     * See https://github.com/broadinstitute/gatk/issues/4724.
     */
    public static String getCanonicalPath(final String filename) {
        Utils.nonEmpty(filename);
        return getCanonicalPath(new File(filename));
    }
}
