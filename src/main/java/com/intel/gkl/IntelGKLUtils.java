/*
 * The MIT License
 *
 * Copyright (c) 2016 Intel Corporation
 *
 * 	Permission is hereby granted, free of charge, to any person
 * 	obtaining a copy of this software and associated documentation
 * 	files (the "Software"), to deal in the Software without
 * 	restriction, including without limitation the rights to use,
 * 	copy, modify, merge, publish, distribute, sublicense, and/or
 * 	sell copies of the Software, and to permit persons to whom the
 * 	Software is furnished to do so, subject to the following
 * 	conditions:
 *
 * 	The above copyright notice and this permission notice shall be
 * 	included in all copies or substantial portions of the
 * 	Software.
 *
 * 	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY
 * 	KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
 * 	WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * 	PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
 * 	COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * 	LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
 * 	OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * 	SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

package com.intel.gkl;

import com.intel.gkl.NativeLibraryLoader;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.gatk.nativebindings.NativeLibrary;

import java.io.File;

/**
 * Provides utilities used by the GKL library.
 */
public final class IntelGKLUtils implements NativeLibrary {
    private final static Logger logger = LogManager.getLogger(com.intel.gkl.IntelGKLUtils.class);
    private static final String NATIVE_LIBRARY_NAME = "gkl_utils";
    private static boolean initialized = false;

    /**
     * Loads the native library, if it is supported on this platform. <p>
     * Returns false if the native library cannot be loaded for any reason. <br>
     *
     * @param tempDir  directory where the native library is extracted or null to use the system temp directory
     * @return  true if the native library is supported and loaded, false otherwise
     */
    @Override


    public synchronized boolean load(File tempDir) {

        if (!NativeLibraryLoader.load(tempDir, NATIVE_LIBRARY_NAME)) {
            System.out.println("NATIVE LIBARARY LOADER FAILED LOADING TO "+tempDir);
            System.err.println("NATIVE LIBARARY LOADER FAILED LOADING TO "+tempDir);
            logger.error("NATIVE LIBARARY LOADER FAILED LOADING TO "+tempDir);
            return false;
        }
        if (!initialized) {
            initialized = true;
        }
        return true;
    }

    private static final String TEST_RESOURCES_PATH = System.getProperty("user.dir") + "/src/test/resources/";
    private static final String TEST_RESOURCES_ABSPATH = new File(TEST_RESOURCES_PATH).getAbsolutePath() + "/";

    public static String pathToTestResource(String filename) {
        return TEST_RESOURCES_ABSPATH + filename;
    }


    public boolean getFlushToZero() {
        return getFlushToZeroNative();
    }

    public void setFlushToZero(boolean value) {
        setFlushToZeroNative(value);
    }

    public boolean isAvxSupported() {
        return isAvxSupportedNative();
    }

    public boolean isAvx2Supported() {
        return isAvx2SupportedNative();
    }

    public boolean isAvx512Supported() {
        return isAvx512SupportedNative();
    }

    /**
     * Get the number of available OpenMP threads, when enabled.
     * @return  Number of available threads, as reported by OpenMP.
     */
    public int getAvailableOmpThreads() {
        return getAvailableOmpThreadsNative();
    }

    private native boolean getFlushToZeroNative();
    private native void setFlushToZeroNative(boolean value);
    private native boolean isAvxSupportedNative();
    private native boolean isAvx2SupportedNative();
    private native boolean isAvx512SupportedNative();
    private native int getAvailableOmpThreadsNative();
}
