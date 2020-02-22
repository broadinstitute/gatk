package org.broadinstitute.hellbender.utils.pairhmm;

import org.broadinstitute.hellbender.utils.Utils;

import java.io.*;

public abstract class BinaryTableWriter<R> implements AutoCloseable {

    private final DataOutputStream dataOut;

    private String path;

    private long counter;

    public BinaryTableWriter(final File file) throws FileNotFoundException {
        this(new FileOutputStream(file), file.toString());
    }

    public BinaryTableWriter(final OutputStream out) {
        this(out, null);
    }

    public BinaryTableWriter(final OutputStream out, final String path) {
        dataOut = new DataOutputStream(Utils.nonNull(out));
        this.path = path;
    }

    public String getPath() {
        return path;
    }

    protected abstract void writeRecord(final R record, final DataOutput output)
            throws IOException;

    public void write(final R record) throws IOException {
        writeRecord(record, dataOut);
        counter++;
    }

    public void writeAll(final Iterable<? extends R> records) throws IOException {
        for (final R record : records) {
            writeRecord(record, dataOut);
            counter++;
        }
    }

    @Override
    public void close() throws IOException {
        dataOut.close();
    }

    public long getCounter() {
        return counter;
    }
}
