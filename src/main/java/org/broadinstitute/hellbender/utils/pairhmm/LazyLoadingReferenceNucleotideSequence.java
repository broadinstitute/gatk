package org.broadinstitute.hellbender.utils.pairhmm;

import htsjdk.samtools.SAMSequenceRecord;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.util.concurrent.*;

public class LazyLoadingReferenceNucleotideSequence implements NucleotideSequence {
    private final String id;
    private final long length;
    private int bufferSize;
    private byte[] buffer;
    private byte[] previousBuffer;
    private final ReferenceDataSource dataSource;
    private long bufferStart;
    private long bufferEnd;


    public LazyLoadingReferenceNucleotideSequence(final ReferenceDataSource dataSource, final String id, final long length, final int bufferSize) {
        this.dataSource = dataSource;
        this.id = id;
        this.length = length;
        this.bufferSize = bufferSize < 128 ? 128 : bufferSize;
        this.buffer = null;
        this.previousBuffer = null;
        this.bufferStart = -1;
        this.bufferEnd = -1;
    }

    private void loadStartingAt(final long start) {
        final long maximumEnd = start + bufferSize - 1;
        final long newEnd = maximumEnd > length ? length : maximumEnd;
        previousBuffer = buffer != null && bufferStart == start - buffer.length ? buffer : null; // only cache previous if it is contiguous.
        buffer = dataSource.queryAndPrefetch(id, start, newEnd).getBases();
        bufferStart = start;
        bufferEnd = newEnd;
    }

    public static LazyLoadingReferenceNucleotideSequence of(final ReferenceDataSource dataSource, final String contigId, final int bufferSize) {
        Utils.nonNull(dataSource);
        Utils.nonNull(contigId);
        ParamUtils.isPositive(bufferSize, "buffer size");
        final SAMSequenceRecord sequenceRecord = dataSource.getSequenceDictionary().getSequence(contigId);
        if (sequenceRecord == null) {
            throw new IllegalArgumentException("there is not such sequence/contig id " + contigId + " in input data-source");
        }
        return new LazyLoadingReferenceNucleotideSequence(dataSource, contigId, sequenceRecord.getSequenceLength(), bufferSize);
    }

    @Override
    public long length() {
        return length;
    }

    public byte byteAt(final long location) {
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
