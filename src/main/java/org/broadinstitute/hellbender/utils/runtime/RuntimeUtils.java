package org.broadinstitute.hellbender.utils.runtime;

import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.net.URL;
import java.util.jar.Manifest;

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

    /**
     * @param clazz class to use when looking up the Implementation-Title
     * @return The name of this toolkit, uses "Implementation-Title" from the
     *         jar manifest of the given class, or (if that's not available) the package name.
     */
    public static String getToolkitName(Class<?> clazz) {
        final String implementationTitle = clazz.getPackage().getImplementationTitle();
        return implementationTitle != null ? implementationTitle : clazz.getPackage().getName();
    }

    /**
     * @return load the manifest file associated with the given class or null if loading fails
     */
    public static Manifest getManifest(Class<?> clazz) {
        final URL resourceURL = clazz.getResource(clazz.getSimpleName() + ".class");
        if( resourceURL == null) {
            return null;
        }
        final String classPath = resourceURL.toString();
        if (classPath.startsWith("jar")) {
            final String manifestPath = classPath.substring(0, classPath.lastIndexOf("!") + 1) + "/META-INF/MANIFEST.MF";
            try ( final InputStream manifestStream = new URL(manifestPath).openStream() ) {
                return new Manifest(manifestStream);
            } catch (IOException e) {
                return null;
            }
        } else {
            return null;
        }
    }

    /**
     * @return get the implementation version of the given class
     */
    public static String getVersion(Class<?> clazz){
        String versionString = clazz.getPackage().getImplementationVersion();
        return versionString != null ? versionString : "Unavailable";
    }
}
