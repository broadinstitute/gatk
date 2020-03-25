package org.broadinstitute.hellbender.tools.dragstr;

import htsjdk.samtools.SAMSequenceRecord;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

/**
 * Buffered access to long reference sequences.
 * <p>
 * Allows to access the whole contig sequence down to a single base at a time without worrying about loading the next frame.
 * </p>
 * Allows contiguous access of individual bases for down to a single base retrieval loading more bases from the underlying sequence as needed.
 * <p>
 *     It keep a copy of the previous loaded section in case there some additional bases up-stream are requested.
 * </p>
 * <p>
 *     This complements any buffering that might be provided by the underlying {@link ReferenceDataSource}.
 * </p>
 */
final class BufferedReferenceBases {
    private final String contigID;
    private final long length;
    private int bufferSize;
    private byte[] buffer;
    private byte[] previousBuffer;
    private final ReferenceDataSource dataSource;
    private long bufferStart;
    private long bufferEnd;

    private BufferedReferenceBases(final ReferenceDataSource dataSource, final String contigID, final long length, final int bufferSize) {
        this.dataSource = dataSource;
        this.contigID = contigID;
        this.length = length;
        this.bufferSize = (int) Math.min(bufferSize , length);
        this.buffer = null;
        this.previousBuffer = null;
        this.bufferStart = -1;
        this.bufferEnd = -1;
    }

    public String contigID() {
        return contigID;
    }

    private void loadStartingAt(final long start) {
        final long maximumEnd = start + bufferSize - 1;
        final long newEnd = Math.min(maximumEnd, length);
        previousBuffer = buffer != null && bufferStart == start - buffer.length ? buffer : null; // only cache previous if it is contiguous.
        buffer = dataSource.queryAndPrefetch(contigID, start, newEnd).getBases();
        bufferStart = start;
        bufferEnd = newEnd;
    }

    public static BufferedReferenceBases of(final ReferenceDataSource dataSource, final String contigId, final int bufferSize) {
        Utils.nonNull(dataSource);
        Utils.nonNull(contigId);
        ParamUtils.isPositive(bufferSize, "buffer size");
        final SAMSequenceRecord sequenceRecord = dataSource.getSequenceDictionary().getSequence(contigId);
        if (sequenceRecord == null) {
            throw new IllegalArgumentException("there is not such sequence/contig id " + contigId + " in input data-source");
        }
        return new BufferedReferenceBases(dataSource, contigId, sequenceRecord.getSequenceLength(), bufferSize);
    }

    public long length() {
        return length;
    }

    final int copyBytesAt(final long position, final byte[] dest, final int offset, final int requestLength) {
        final long lastPosition = Math.min(length, position + requestLength - 1);
        long p;
        int o;
        for (p = position, o = offset; p <= lastPosition; ++p, ++o) {
            dest[o] = byteAt(p);
        }
        return (int) (lastPosition - position + 1);
    }

    byte byteAt(final long location) {
        if (location < bufferStart) {
            if (location < 0) {
                throw new IllegalArgumentException("invalid location outside range: " + location);
            } else if (previousBuffer != null && location >= bufferStart - previousBuffer.length) {
                return previousBuffer[(int) (location - (bufferStart - previousBuffer.length))];
            } else {
                loadStartingAt(location);
                return buffer[0];
            }
        } else if (location > bufferEnd) {
            if (location > length) {
                throw new IllegalArgumentException("invalid location");
            } else {
                loadStartingAt(location);
                return buffer[0];
            }
        } else {
            return buffer[(int) (location - bufferStart)];
        }
    }

}
