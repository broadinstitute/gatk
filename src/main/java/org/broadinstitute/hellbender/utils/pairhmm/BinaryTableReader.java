package org.broadinstitute.hellbender.utils.pairhmm;

import org.broadinstitute.hellbender.exceptions.UserException;

import javax.swing.*;
import java.io.*;
import java.util.*;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;

public abstract class BinaryTableReader<R> implements AutoCloseable, Iterator<R> {

    private final PushbackDataInput dtInput;

    private final String inputName;

    private final Runnable closeAction;

    protected BinaryTableReader(final InputStream in, final String source) {
        final DataInputStream dataInputStream = new DataInputStream(in);
        dtInput = new ByteBufferPushbackDataInput(dataInputStream, 1024);
        inputName = source;
        this.closeAction = () -> {
            try {
                dataInputStream.close();
            } catch (final IOException ex) {
                throw source != null ? new UserException.CouldNotReadInputFile(source)
                                     : new UserException.CouldNotReadInputFile("unknown source");
            }
        };
    }

    protected BinaryTableReader(final DataInput in) {
        this.dtInput = new ByteBufferPushbackDataInput(in, 1024);
        this.inputName = null;
        this.closeAction = () -> {};
    }

    protected abstract R readRecord(final PushbackDataInput input)
        throws IOException;

    protected boolean hasMoreRecords(final PushbackDataInput input)
        throws IOException {
        return !input.eof();
    }

    public Stream<R> stream() {
        return StreamSupport.stream(Spliterators.spliteratorUnknownSize(new Iterator<R>() {
            private R next;
            @Override
            public boolean hasNext() {
                if (next != null) {
                    return false;
                }
                try {
                    next = readRecord(dtInput);
                    return next != null;
                } catch (final EOFException ex) {
                    return false;
                } catch (final IOException ex) {
                    throw new UncheckedIOException(ex);
                }
            }

            @Override
            public R next() {
                if (next != null || hasNext()) {
                    final R result = next;
                    next = null;
                    return result;
                } else {
                    throw new NoSuchElementException();
                }
            }
        },0), false);
    }

    public void close() throws IOException {
        try {
            closeAction.run();
        } catch (final Exception ex) {
            if (ex.getCause() instanceof IOException) {
                throw (IOException) ex.getCause();
            } else {
                throw ex;
            }
        }
    }

    public List<R> readAll() throws IOException {
        final List<R> result = new ArrayList<>();
        while (hasMoreRecords(dtInput)) {
            final R record = readRecord(dtInput);
            if (record == null) {
                return result;
            }
            result.add(record);
        }
        return result;
    }

    public boolean hasNext() {
        try {
            return hasMoreRecords(dtInput);
        } catch (final IOException ex) {
            throw new UncheckedIOException(ex);
        }
    }

    public R next() {
        try {
            if (hasMoreRecords(dtInput)) {
                return readRecord(dtInput);
            } else {
                throw new NoSuchElementException();
            }
        } catch (final IOException ex) {
            throw new UncheckedIOException(ex);
        }
    }
}
