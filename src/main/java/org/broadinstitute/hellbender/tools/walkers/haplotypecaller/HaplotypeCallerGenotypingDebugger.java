package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.io.IOException;
import java.io.PrintStream;
import java.nio.file.Files;

/**
 * A short helper class that manages a singleton debug stream for HaplotypeCaller genotyping information that is useful for debugging.
 *
 * In order to use simply call initialize() providing a location for an output path, then call isEnabled() to evaluate if the processing
 * for printing debug statements should be performed, and call println() to output to the file in a synchronized fashion.
 */
public class HaplotypeCallerGenotypingDebugger{
    private static PrintStream genotyperDebugOutStream;

    public static void initialize(final String debugLocation) {
        try {
            genotyperDebugOutStream = new PrintStream(Files.newOutputStream(IOUtils.getPath(debugLocation)));
        } catch (IOException e) {
            throw new UserException.CouldNotCreateOutputFile(debugLocation, "Provided argument for genotyper debug location could not be created");
        }
    }

    // Is the debugger enabled
    public static boolean isEnabled() {
        return genotyperDebugOutStream != null;
    }

    // Print the provided text to the debugger if it exists
    public static synchronized void println(final String debug) {
        if (genotyperDebugOutStream != null) {
            genotyperDebugOutStream.println(debug);
        }
    }

    // Closes out the debug output stream if necessary
    public static void close() {
        if (genotyperDebugOutStream != null) {
            genotyperDebugOutStream.close();
        }
    }
}
