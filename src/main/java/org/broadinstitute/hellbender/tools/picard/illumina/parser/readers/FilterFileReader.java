package org.broadinstitute.hellbender.tools.picard.illumina.parser.readers;

import org.broadinstitute.hellbender.utils.UnsignedTypeUtil;

import java.io.File;
import java.nio.ByteBuffer;
import java.util.*;

import static java.lang.Integer.toHexString;
import static org.broadinstitute.hellbender.tools.picard.illumina.parser.readers.MMapBackedIteratorFactory.getByteIterator;
import static org.broadinstitute.hellbender.utils.UnsignedTypeUtil.uIntToLong;

/**
 * Illumina uses an algorithm described in "Theory of RTA" that determines whether or not a cluster passes filter("PF") or not.
 * These values are written as sequential bytes to Filter Files.  The structure of a filter file is as follows:
 * Bytes 0-3  : 0
 * Bytes 4-7  : unsigned int version
 * Bytes 8-11 : unsigned int numClusters
 */
public class FilterFileReader implements Iterator<Boolean> {
    /**
     * Number of bytes in the files header that will be skipped by the iterator
     */
    private static final int HEADER_SIZE = 12;

    /**
     * Expected Version
     */
    public final int EXPECTED_VERSION = 3;

    /**
     * Iterator over each cluster in the FilterFile
     */
    private final BinaryFileIterator<Byte> bbIterator;

    /**
     * Version number found in the FilterFile, this should equal 3
     */
    public final int version;

    /**
     * The number of cluster's pf values stored in this file
     */
    public final long numClusters;

    /**
     * Byte representing a cluster failing filter(not a PF read), we test this exactly at
     * the moment but technically the standard  may be to check only lowest significant bit
     */
    private final static byte FailedFilter = 0x00;

    /**
     * Byte representing a cluster passing filter(a PF read), we test this exactly at
     * the moment but technically the standard  may be to check only lowest significant bit
     */
    private final static byte PassedFilter = 0x01;

    /**
     * The index of the current cluster within the file
     */
    private int currentCluster;

    public FilterFileReader(final File file) {
        bbIterator = getByteIterator(HEADER_SIZE, file);
        final ByteBuffer headerBuf = bbIterator.getHeaderBytes();

        for (int i = 0; i < 4; i++) {
            final byte b = headerBuf.get();
            if (b != 0) {
                throw new IlluminaReaderException("The first four bytes of a Filter File should be 0 but byte " + i + " was " + b + " in file " + file.getAbsolutePath());
            }
        }

        version = headerBuf.getInt();
        if (version != EXPECTED_VERSION) {
            throw new IlluminaReaderException("Expected version is " + EXPECTED_VERSION + " but version found was " + version + " in file " + file.getAbsolutePath());
        }

        numClusters = uIntToLong(headerBuf.getInt());
        bbIterator.assertTotalElementsEqual(numClusters);

        currentCluster = 0;
    }

    public boolean hasNext() {
        return currentCluster < numClusters;
    }

    public Boolean next() {
        final byte value = bbIterator.next();
        currentCluster += 1;
        if (value == PassedFilter) {
            return true;
        } else if (value == FailedFilter) {
            return false;
        } else {
            String hexVal = toHexString(value);
            hexVal = (hexVal.length() < 2 ? "0x0" : "0x") + hexVal;
            throw new IlluminaReaderException("Didn't recognized PF Byte (" + hexVal + ")" + " for element (" + currentCluster + ") in file(" + bbIterator.getFile().getAbsolutePath() + ")");
        }
    }

    public void skipRecords(final int numToSkip) {
        bbIterator.skipElements(numToSkip);
    }

    public void remove() {
        throw new UnsupportedOperationException();
    }
}
