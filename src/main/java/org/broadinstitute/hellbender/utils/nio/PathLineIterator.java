package org.broadinstitute.hellbender.utils.nio;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.IOException;
import java.nio.charset.CharacterCodingException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Iterator;
import java.util.Spliterator;
import java.util.function.Consumer;
import java.util.stream.Stream;

/**
 * Iterate through the lines of a Path. Works for anything you can point
 * a Path to, such as a normal file or anything you have an NIO provider
 * for (e.g. GCS).
 */
public class PathLineIterator implements AutoCloseable, Iterable<String> {
    private final Stream<String> lines;

    /**
     * Returns an iterator so you can iterate over the lines in the text file like so:
     * for (String line: new PathLineIterator(path)) {
     *   // do something with the line
     * }
     *
     * It's also closeable so you can close it when done, or use it in a try-with-resources
     * to close it automatically.
     *
     * @param path path to a text file.
     */
    public PathLineIterator(final Path path) {
        try {
            lines = Files.lines(Utils.nonNull(path, "path shouldn't be null"));
        }
        catch (final CharacterCodingException ex ) {
            throw new UserException("Error detected in file character encoding.  Possible inconsistent character encodings within the file: " + path.toUri().toString(), ex);
        }
        catch (final IOException x) {
            throw new UserException("Error reading " + path.toUri().toString(), x);
        }
    }

    @Override
    public void close() {
        lines.close();
    }

    @Override
    public Iterator<String> iterator() {
        return lines.iterator();
    }

    @Override
    public void forEach(Consumer<? super String> action) {
        lines.forEach(action);
    }

    @Override
    public Spliterator<String> spliterator() {
        return lines.spliterator();
    }
}
