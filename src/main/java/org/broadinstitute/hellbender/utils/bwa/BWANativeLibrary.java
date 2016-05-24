package org.broadinstitute.hellbender.utils.bwa;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.NativeUtils;

/**
 * A class to help with loading the native JBWA code. All clients who want to use JBWA need to call
 * BWANativeLibrary.loadNativeJBWALibrary() first.
 */
public class BWANativeLibrary {
    public static void load() {
        // Load native library in each task VM
        final String libraryPath = NativeUtils.runningOnLinux() ? "/lib/libbwajni.so" : "/lib/libbwajni_mac.jnilib";
        if (! NativeUtils.loadLibraryFromClasspath(libraryPath)){
            throw new UserException( "jbwa library was not loaded. This could be due to a configuration error, or your system might not support it.");
        }
    }
}
