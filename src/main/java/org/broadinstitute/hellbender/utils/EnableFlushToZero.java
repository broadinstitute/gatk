package org.broadinstitute.hellbender.utils;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.io.Resource;

import java.io.*;

/**
 * Class to enable Flush To Zero (FTZ) mode for floating-point operations on x86 processors
 *
 * Enabling FTZ provides better performance when handling floating-point underflow conditions.
 * When FTZ is enabled, floating-point SIMD instructions (SSE/AVX) that result in a underflow will return a zero,
 * instead of a denormalized floating-point number.
 * Floating-point SIMD instructions are used in the native AVX PairHMM implementation as well as code generated
 * by the JVM.
 *
 */
public final class EnableFlushToZero {
    private static final Logger logger = LogManager.getLogger(EnableFlushToZero.class);
    private static final Boolean runningOnMac = System.getProperty("os.name", "unknown").toLowerCase().startsWith("mac");
    private static final String path = "/lib/libEnableFlushToZero" + (runningOnMac ? ".dylib" : ".so");

    public static void enable() {
        try {
            // copy the shared library from the jar file to a temp file
            File temp = IOUtils.writeTempResource(new Resource(path, EnableFlushToZero.class));
            // delete the temp file when the JVM exits
            temp.deleteOnExit();
            // load the shared library, which enables FTZ
            System.load(temp.getAbsolutePath());
            logger.info("Enabled Flush To Zero mode.");
        } catch (IllegalArgumentException|UnsatisfiedLinkError e) {
            logger.warn("Cannot load '" + path + "' from jar file.");
            logger.warn("Cannot enable Flush To Zero mode. RESULTS MAY BE DIFFERENT THAN EXPECTED compared to Java with Flush To Zero enabled.");
        }
    }
}
