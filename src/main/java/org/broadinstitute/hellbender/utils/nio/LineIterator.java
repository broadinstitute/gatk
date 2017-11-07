package org.broadinstitute.hellbender.utils.nio;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Iterator;
import java.util.Spliterator;
import java.util.function.Consumer;
import java.util.stream.Stream;
import org.broadinstitute.hellbender.exceptions.UserException;

/**
 * Iterate through the lines of a file.
 */
public class LineIterator implements AutoCloseable, Iterable<String> {
  private final Stream<String> lines;

  /**
   * Returns an iterator so you can iterate over the lines in the text file like so:
   * for (String line: new Utils.LineIterator(path)) {
   *   // do something with the line
   * }
   *
   * It's also closeable so you can close it when done, or use it in a try-with-resources
   * to close it automatically.
   *
   * @param path path to a text file.
   * @throws UserException if we cannot open the file for reading.
   */
  public LineIterator(Path path) throws UserException {
    try {
      lines = Files.lines(path);
    } catch (IOException x) {
      throw new UserException("Error reading " + path.toString() + ": " + x.getMessage());
    }
  }

  /**
   * @InheritDoc
   */
  @Override
  public void close() {
    lines.close();
  }

  /**
   * @return an Iterator over the lines (as Strings).
   */
  @Override
  public Iterator<String> iterator() {
    return lines.iterator();
  }

  /**
   * @InheritDoc
   */
  @Override
  public void forEach(Consumer<? super String> action) {
    lines.forEach(action);
  }

  /**
   * @InheritDoc
   */
  @Override
  public Spliterator<String> spliterator() {
    return lines.spliterator();
  }
}
