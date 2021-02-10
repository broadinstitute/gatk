package org.broadinstitute.hellbender.cmdline.argumentcollections;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.exceptions.UserException;

import java.io.IOException;
import java.io.Serializable;
import java.nio.file.Files;
import java.nio.file.StandardOpenOption;

/**
 * To Use this class add a @ArgumentCollection variable like so:
 *
 * <code>
 *
 * @ArgumentCollection public final out = new SimpleOutput();
 * </code>
 * <p>
 * and in the method <code>onTraversalSuccess</code> instead of just outputting to the terminal, also call
 *
 * <code>
 * out.writeToOutput(value)
 * </code>
 * <p>
 * where <code>value</code> is the object to be written to the file.
 * <p>
 * The code will only attempt to write to the output if -o (or --output) has been assigned on the commandline.
 */
public class OptionalTextOutputArgumentCollection implements Serializable {
    private static final long serialVersionUID = 1L;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "Optional file for the output value.",
            optional = true)
    public GATKPath output = null;

    /**
     * Print value (value.toString()) to the output path, if output is not null
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
     * Print value (value.toString() + "\n") to the output path, if output is not null
     */
    public void println(Object value) {
        print(value);
        print("\n");
    }
}
