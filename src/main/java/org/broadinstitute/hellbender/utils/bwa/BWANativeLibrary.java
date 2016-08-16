package org.broadinstitute.hellbender.utils.bwa;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.NativeUtils;

/**
 * A class to help with loading the native JBWA code. All clients who want to use JBWA need to call
 * BWANativeLibrary.loadNativeJBWALibrary() first.
 */
public class BWANativeLibrary {
    /**
     * load the native JBWA Library, this should be called before using any methods from {@link com.github.lindenb.jbwa.jni}
     * @throws org.broadinstitute.hellbender.exceptions.UserException.HardwareFeatureException if the library could not be loaded
     */
    public static void load() {
        // Load native library in each task VM
        final String libraryPath = getLibraryNameForSystem();
        if (!NativeUtils.loadLibraryFromClasspath(libraryPath)) {
            throw new UserException.HardwareFeatureException("jbwa library was not loaded. This could be due to a configuration error, or your system might not support it.");
        }
    }

    private static String getLibraryNameForSystem() {
        if (NativeUtils.runningOnPPCArchitecture()) {
            return "/libbwajni_ppc64.so";
        } else if (NativeUtils.runningOnLinux()) {
            return "/libbwajni.so";
        } else {
            return "/libbwajni.jnilib";
        }
    }
}
