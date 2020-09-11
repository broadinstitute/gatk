package org.broadinstitute.hellbender.utils.utils;

import org.broadinstitute.hellbender.utils.Utils;

import java.io.*;

/**
 * Abstract base class to write into binary formatted table files.
 * <p>
 *     You need to provide a procedure to encode the record from the data input stream overriding
 *     {@link #writeRecord(Object, DataOutput) writeRecord},
 *     and you are good to go.
 * </p>
 * @param <R> the record/row type.
 */
public abstract class BinaryTableWriter<R> implements AutoCloseable {

    private final DataOutputStream dataOut;

    private String path;

    private long counter;

    public BinaryTableWriter(final OutputStream out, final String path) {
        dataOut = new DataOutputStream(Utils.nonNull(out));
        this.path = path;
    }

    public String getPath() {
        return path;
    }

    protected abstract void writeRecord(final R record, final DataOutput output)
            throws IOException;

    public final void write(final R record) throws IOException {
        writeRecord(record, dataOut);
        counter++;
    }

    public final void writeAll(final Iterable<? extends R> records) throws IOException {
        for (final R record : records) {
            writeRecord(record, dataOut);
            counter++;
        }
    }

    @Override
    public final void close() throws IOException {
        dataOut.close();
    }

    public long getCounter() {
        return counter;
    }
}
