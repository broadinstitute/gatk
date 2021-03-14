package org.broadinstitute.hellbender.utils;

import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.io.*;

/**
 * Abstract file writing class for record tables stored in binary format.
 * @param <R> record type.
 */
public abstract class BinaryTableWriter<R> implements AutoCloseable {

    protected final DataOutputStream dataOut;

    private final ByteCounterOutputStream byteCounter;

    private String path;

    private long counter;

    protected BinaryTableWriter(final OutputStream out, final String path) {
        Utils.nonNull(out);
        byteCounter = new ByteCounterOutputStream(out);
        dataOut = new DataOutputStream(byteCounter);
        this.path = path;
    }

    public long offset() {
        return byteCounter.count();
    }

    public String getPath() {
        return path;
    }

    protected abstract void writeRecord(final R record, final DataOutput output)
            throws IOException;

    public void write(final R record) throws IOException {
        writeRecord(Utils.nonNull(record), dataOut);
        counter++;
    }

    public void writeAll(final Iterable<? extends R> records) throws IOException {
        for (final R record : Utils.nonNull(records)) {
            writeRecord(Utils.nonNull(record), dataOut);
            counter++;
        }
    }

    @Override
    public void close() throws IOException {
        dataOut.close();
    }

    public void flush() throws IOException {
        dataOut.flush();
    }

    public long numberOfRecords() {
        return counter;
    }

    private static class ByteCounterOutputStream extends OutputStream {
        private final OutputStream out;
        private long count;

        long count() {
            return count;
        }

        private ByteCounterOutputStream(OutputStream out) {
            this.out = Utils.nonNull(out);
            this.count = 0;
        }

        @Override
        public void write(int b) throws IOException {
            Utils.nonNull(b);
            out.write(b);
            count++;
        }

        @Override
        public void write(byte[] b) throws IOException {
            Utils.nonNull(b);
            out.write(b);
            count += b.length;
        }

        @Override
        public void write(byte[] b, int off, int len) throws IOException {
            Utils.nonNull(b);
            ParamUtils.isValidArrayOffsetAndRangeLength(off, len, b.length, "input buffer");
            out.write(b, off, len);
            count += len;
        }

        @Override
        public void flush() throws IOException {
            out.flush();
        }

        @Override
        public void close() throws IOException {
            out.close();
        }
    }
}
