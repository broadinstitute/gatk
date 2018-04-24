package org.broadinstitute.hellbender.utils;

import org.apache.commons.lang3.SystemUtils;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.io.Resource;

import java.io.File;

/**
 * Utilities to provide architecture-dependent native functions
 */
public final class NativeUtils {

    private NativeUtils(){}

    /**
     * boolean : returns whether the current platform is PowerPC, LE only
     */
    public static boolean runningOnPPCArchitecture() {
        return runningOnArchitecture("ppc64le");
    }

    /**
     * @param architecture Architecture to look for
     * @return true if we're running on the specified architecture, otherwise false
     */
    public static boolean runningOnArchitecture( final String architecture ) {
        Utils.nonNull(architecture);
        final String archtectureString = SystemUtils.OS_ARCH;

        return archtectureString != null && archtectureString.equals(architecture);
    }

    /**
     * @return true if we're running on a Mac operating system, otherwise false
     */
    public static boolean runningOnMac() {
        return SystemUtils.IS_OS_MAC;
    }

    /**
     * @return true if we're running on a Linux operating system, otherwise false
     */
    public static boolean runningOnLinux() {
        return SystemUtils.IS_OS_LINUX;
    }

    /**
     * Loads a native library stored on our classpath by extracting it to a temporary location and invoking
     * {@link System#load} on it
     *
     * @param libraryPathInJar absolute path to the library file on the classpath
     * @return true if the library was extracted and loaded successfully, otherwise false
     */
    public static boolean loadLibraryFromClasspath( final String libraryPathInJar ) {
        Utils.nonNull(libraryPathInJar);
        Utils.validateArg(libraryPathInJar.startsWith("/"), "library path in jar must be absolute");

        try {
            final File extractedLibrary = IOUtils.writeTempResourceFromPath(libraryPathInJar, NativeUtils.class);
            extractedLibrary.deleteOnExit();
            System.load(extractedLibrary.getAbsolutePath());
        }
        catch ( UnsatisfiedLinkError | UserException e ) {
            return false;
        }

        return true;
    }
}
