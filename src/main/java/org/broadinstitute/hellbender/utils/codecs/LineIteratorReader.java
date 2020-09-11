package org.broadinstitute.hellbender.utils.codecs;

import htsjdk.tribble.readers.LineIterator;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.IOException;
import java.io.Reader;

/**
 * A reader wrapper around a LineIterator.
 *
 * <p>
 *     This class is useful to implement <i>codecs</i> for formats for which a {@link Reader} based parser already
 *     exists.
 * </p>
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public class LineIteratorReader extends Reader {

    private final LineIterator lineIterator;

    private final StringBuffer buffer;
    private int nextChar = 0;

    public LineIteratorReader(final LineIterator lineIterator) {
        this.lineIterator = Utils.nonNull(lineIterator);
        buffer = new StringBuffer();
        nextChar = 0;
    }

    /**
     * Reads lines from the buffer, and if empty from the underlying line-iterator.
     *
     * @param cbuf the destination buffer.
     * @param off  the start position in the buffer to copy characters into.
     * @param len  number of characters to read.
     * @return never less than 0 as this buffer never ends.
     */
    @Override
    public int read(final char[] cbuf, final int off, final int len)  {

        if (len == 0) {
            return 0;
        }
        int readSoFar = 0;
        while (readSoFar < len) {
            final int inBufferCharCount = buffer.length() - nextChar;
            if (inBufferCharCount == 0) {
                if (lineIterator.hasNext()) {
                    buffer.setLength(0);
                    nextChar = 0;
                    buffer.append(lineIterator.next()).append('\n');
                } else {
                    break;
                }
            }
            final int inBuffer = Math.min(len - readSoFar, buffer.length() - nextChar);
            for (int i = 0; i < inBuffer; ++i) {
                cbuf[off + (readSoFar++)] = buffer.charAt(nextChar++);
            }
        }
        return readSoFar == 0 ? -1 : readSoFar;
    }

    @Override
    public void close() throws IOException {
        // nothing to do here.
    }
}
