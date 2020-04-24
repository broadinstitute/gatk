package org.broadinstitute.hellbender.utils.pairhmm;

import java.io.*;

public abstract class BinaryTableWriter<R> implements AutoCloseable {

    protected final DataOutputStream dataOut;

    private final ByteCounterOutputStream byteCounter;

    private String path;

    private long counter;

    public BinaryTableWriter(final File file) throws FileNotFoundException {
        this(new FileOutputStream(file), file.toString());
    }

    public BinaryTableWriter(final OutputStream out) {
        this(out, null);
    }

    public BinaryTableWriter(final OutputStream out, final String path) {
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

    private static class ByteCounterOutputStream extends OutputStream {
        private final OutputStream out;
        private long count;

        long count() {
            return count;
        }

        private ByteCounterOutputStream(OutputStream out) {
            this.out = out;
            this.count = 0;
        }

        @Override
        public void write(int b) throws IOException {
            out.write(b);
            count++;
        }

        @Override
        public void write(byte[] b) throws IOException {
            out.write(b);
            count += b.length;
        }

        @Override
        public void write(byte[] b, int off, int len) throws IOException {
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
