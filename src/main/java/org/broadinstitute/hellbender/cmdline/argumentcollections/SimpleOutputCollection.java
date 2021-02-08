package org.broadinstitute.hellbender.cmdline.argumentcollections;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.engine.GATKPath;

import java.io.IOException;
import java.io.Serializable;
import java.nio.file.Files;
import java.nio.file.StandardOpenOption;

/**
 * To Use this class add a @ArgumentCollection variable like so:
 *
 * <code>
 * @ArgumentCollection
 * public final out = new SimpleOutput();
 * </code>
 *
 * and in the method <code>onTraversalSuccess</code> instead of just outputting to the terminal, also call
 *
 * <code>
 *     out.writeToOutput(value)
 * </code>
 *
 * where <code>value</code> is the object to be written to the file.
 *
 * The code will only attempt to write to the output if -o (or --output) has been assigned on the commandline.
 *
 */
public class SimpleOutputCollection implements Serializable {
    private static final long serialVersionUID = 1L;

    @Argument(shortName = "o", doc = "Optional output file for the output value.", optional = true)
    public GATKPath output = null;


    public void writeToOutput(Object value)  {
        if (output != null) {
            try {
                Files.write(output.toPath(), value.toString().getBytes());
                Files.write(output.toPath(), "\n".getBytes(), StandardOpenOption.APPEND);
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
    }
}
