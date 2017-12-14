package org.broadinstitute.hellbender.utils.runtime;

import org.apache.commons.lang3.StringUtils;

import java.io.File;

import org.broadinstitute.hellbender.utils.Utils;

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

    /**
     * Given a Class that is a CommandLineProgram, either GATK or Picard, return a display name suitable for
     * presentation to the user that distinguishes GATK tools from Picard tools by including a " (Picard)" suffix;
     * @param toolClass A CommandLineProgram class object may not be null
     * @return tool display name
     */
    public static String toolDisplayName(final Class<?> toolClass) {
        Utils.nonNull(toolClass, "A valid class is required to get a display name");

        final String picardToolSuffix = " (Picard)";
        final Class<?> picardCommandLineProgramClass = picard.cmdline.CommandLineProgram.class;
        return picardCommandLineProgramClass.isAssignableFrom(toolClass) ?
                toolClass.getSimpleName() + picardToolSuffix :
                toolClass.getSimpleName();
    }

}
