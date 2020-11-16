package org.broadinstitute.hellbender.utils.io;

import htsjdk.tribble.Feature;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.IOException;
import java.io.PrintStream;
import java.util.function.Function;

/**
 * Encodes Tribble {@link Feature}s to a text file, given a feature encoding function
 */
public final class UncompressedFeatureOutputStream<F extends Feature> implements FeatureOutputStream<F> {

    private static final String NEWLINE_CHARACTER = "\n";

    private final PrintStream outputStream;
    private final Function<F, String> encoder;

    /**
     * @param file file to write to
     * @param encoder encoding function mapping a feature to its text record (whithout a newline character)
     */
    public UncompressedFeatureOutputStream(final GATKPath file, final Function<F, String> encoder) {
        Utils.nonNull(file);
        Utils.nonNull(encoder);
        outputStream = new PrintStream(file.getOutputStream());
        this.encoder = encoder;
    }

    /**
     * Writes header to file. Should be called before adding features
     *
     * @param header header text (without final newline character), cannot be null
     */
    public void writeHeader(final String header) {
        Utils.nonNull(header);
        try {
            outputStream.write((header + NEWLINE_CHARACTER).getBytes());
        } catch (final IOException e) {
            throw new GATKException("Error writing header", e);
        }
    }

    /**
     * Adds new feature and writes to file
     *
     * @param feature feature to write, cannot be null
     */
    public void add(final F feature) {
        Utils.nonNull(feature);
        outputStream.println(encoder.apply(feature));
    }

    /**
     * Closes the underlying stream
     */
    @Override
    public void close() {
        outputStream.close();
    }
}
