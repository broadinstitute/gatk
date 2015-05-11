package org.broadinstitute.hellbender.utils.runtime;

import org.apache.commons.lang3.StringUtils;

import java.io.File;

public final class RuntimeUtils {
    public static final String[] PATHS;

    static {
        String path = System.getenv("PATH");
        if (path == null)
            path = System.getenv("path");
        if (path == null) {
            PATHS = new String[0];
        } else {
            PATHS = StringUtils.split(path, File.pathSeparatorChar);
        }
    }

    /**
     * Returns the path to an executable or null if it doesn't exist.
     * @param executable Relative path
     * @return The absolute file path.
     */
    public static File which(String executable) {
        for (String path: PATHS) {
            File file = new File(path, executable);
            if (file.exists())
                return file.getAbsoluteFile();
        }
        return null;
    }
}
