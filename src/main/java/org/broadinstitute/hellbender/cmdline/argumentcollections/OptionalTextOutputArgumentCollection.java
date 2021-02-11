package org.broadinstitute.hellbender.cmdline.argumentcollections;

import com.google.common.annotations.VisibleForTesting;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.exceptions.UserException;

import java.io.IOException;
import java.io.Serializable;
import java.nio.file.Files;
import java.nio.file.StandardOpenOption;

/**
 * An ArgumentCollection with an optional output argument, and utility methods for printing String output to it
 *
 * To use this class add an @ArgumentCollection variable to your tool like so:
 *
 * <code>
 * @ArgumentCollection
 * public final out = new OptionalTextOutputArgumentCollection();
 * </code>
 * <p>
 * and in the method <code>onTraversalSuccess</code> instead of just outputting to the terminal, also call
 *
 * <code>
 * out.println(value)
 * </code>
 * or
 * <code>
 * out.print(value)
 * </code>
 * <p>
 * where <code>value</code> is the object to be written to the file.
 * <p>
 * The code will only attempt to write to the output if -O (or --output) has been specified on the commandline.
 */
public class OptionalTextOutputArgumentCollection implements Serializable {
    private static final long serialVersionUID = 1L;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "Optional output file",
            optional = true)
    @VisibleForTesting
    GATKPath output = null;

    /**
     * @return The -O/--output Path specified on the command line, or null if there was none
     */
    public GATKPath getOutputPath() {
        return output;
    }
    
    /**
     * Prints value (value.toString()) to the output path, if output is not null
     * Overwrites any pre-existing output file rather than appending.
     */
    public void print(Object value) {
        if (output != null) {
            try {
                Files.write(output.toPath(), value.toString().getBytes());
            } catch (IOException e) {
                throw new UserException.CouldNotCreateOutputFile(output.toString(), e);
            }
        }
    }

    /**
     * Prints value (value.toString() + "\n") to the output path, if output is not null
     * Overwrites any pre-existing output file rather than appending.
     */
    public void println(Object value) {
        print(value.toString() + "\n");
    }
}
