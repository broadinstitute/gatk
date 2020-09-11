package org.broadinstitute.hellbender.utils;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;

import javax.swing.*;
import java.io.*;
import java.util.*;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;

/**
 * Abstract base class to read for binary formatted table files.
 * <p>
 *     You need to provide a procedure to decode the record from the data input stream overriding {@link #readRecord(DataInput)},
 *     and you are good to go.
 * </p>
 * @param <R> the record/row type.
 */
public abstract class BinaryTableReader<R> implements AutoCloseable, Iterator<R> {

    private final DataInput dtInput;

    private final Runnable closeAction;

    private R nextRecord;

    public BinaryTableReader(final DataInput in) {
        this.dtInput = Utils.nonNull(in);
        this.closeAction = () -> {};
        try {
            nextRecord = readRecord(dtInput);
        } catch (final EOFException ex) {
            nextRecord = null;
        } catch (final IOException ex) {
            throw new UncheckedIOException(ex);
        }
    }

    protected abstract R readRecord(final DataInput input)
        throws IOException;

    @Override
    public final void close() {
        closeAction.run();
    }

    public final List<R> readAll() {
        final List<R> result = new ArrayList<>();
        while (nextRecord != null) {
            result.add(nextRecord);
            try {
                nextRecord = readRecord(dtInput);
            } catch (final EOFException ex) {
                nextRecord = null;
            } catch (final IOException ex) {
                throw new UncheckedIOException(ex);
            }
        }
        return result;
    }

    public final boolean hasNext() {
        return nextRecord != null;
    }

    public final R next() {
        if (nextRecord != null) {
            final R result = nextRecord;
            try {
                nextRecord = readRecord(dtInput);
            } catch (final EOFException ex) {
                nextRecord = null;
            } catch (final IOException ex) {
                throw new UncheckedIOException(ex);
            }
            return result;
        } else {
            throw new NoSuchElementException();
        }
    }
}
