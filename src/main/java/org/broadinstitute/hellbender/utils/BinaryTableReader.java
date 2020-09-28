package org.broadinstitute.hellbender.utils;

import org.apache.commons.io.input.NullInputStream;
import org.broadinstitute.hellbender.exceptions.UserException;

import java.io.*;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;

/**
 * Abstract base class for readers of table with records stored in binary.
 * @param <R> record type.
 */
public abstract class BinaryTableReader<R> implements AutoCloseable, Iterator<R> {

    private final DataInput dtInput;
    private R next;
    private final Runnable closeAction;

    protected BinaryTableReader(final InputStream in, final String source) {
        final DataInputStream dataInputStream = new DataInputStream(in);
        dtInput = new DataInputStream(in);
        this.closeAction = () -> {
            try {
                dataInputStream.close();
            } catch (final IOException ex) {
                throw source != null ? new UserException.CouldNotReadInputFile(source)
                                     : new UserException.CouldNotReadInputFile("unknown source");
            }
        };
        next = readNextRecord();
    }

    public static <E> BinaryTableReader<E> emptyReader() {
        return new BinaryTableReader<E>(new NullInputStream(0), "null") {
            @Override
            protected E readRecord(final DataInput input) throws EOFException {
                throw new EOFException("reached the end");
            }
        };
    }

    private R readNextRecord() {
        try {
            return readRecord(dtInput);
        } catch (final EOFException ex) {
            return null;
        } catch (final IOException ex) {
            throw new UncheckedIOException(ex);
        }
    }

    protected abstract R readRecord(final DataInput input)
        throws IOException;


    public final Stream<R> stream() {
        return StreamSupport.stream(Spliterators.spliteratorUnknownSize(new Iterator<R>() {
            @Override
            public boolean hasNext() {
                return next != null;
            }

            @Override
            public R next() {
                final R result = next;
                next = readNextRecord();
                return result;
            }
        },0), false);
    }

    public final void close() throws IOException {
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

    public final List<R> readAll() {
        return stream().collect(Collectors.toList());
    }

    public final boolean hasNext() {
        return next != null;
    }

    public final R next() {
        if (next == null) {
            throw new NoSuchElementException();
        } else {
            final R result = next;
            next = readNextRecord();
            return result;
        }
    }
}
