package org.broadinstitute.hellbender.tools.picard.illumina.parser.readers;

import org.broadinstitute.hellbender.utils.UnsignedTypeUtil;

import java.io.File;
import java.nio.ByteBuffer;

import static org.broadinstitute.hellbender.tools.picard.illumina.parser.readers.MMapBackedIteratorFactory.getFloatIterator;
import static org.broadinstitute.hellbender.utils.UnsignedTypeUtil.uIntToLong;

/**
 * The locs file format is one 3 Illumina formats(pos, locs, and clocs) that stores position data exclusively.
 * locs files store position data for successive clusters in 4 byte float pairs, described as follows:
 * bytes 1-4    : (int?) Version number (1)
 * bytes 5-8    : 4 byte float equaling 1.0
 * bytes 9-12   : unsigned int numClusters
 * bytes 13-16: : X coordinate of first cluster (32-bit float)
 * bytes 17-20: : Y coordinate of first cluster (32-bit float)
 * <p>
 * The remaining bytes of the file store the X and Y coordinates of the remaining clusters.
 */

public final class LocsFileReader extends AbstractIlluminaPositionFileReader {
    /**
     * Size of the opening file header, this is skipped by the iterator below
     */
    private static final int HEADER_SIZE = 12;

    /**
     * The first four bytes of a locs file should equal a little endian 1
     */
    private static final int BYTES_1_TO_4 = 1;

    /**
     * The expected version of locs files
     */
    private static final float VERSION = 1.0f;

    /**
     * An iterator over all of the coordinate values in the file, remember next needs to be called
     * twice per coordinate pair
     */
    private BinaryFileIterator<Float> bbIterator;

    /**
     * Total clusters in the file as read in the file header
     */
    private long numClusters;

    /**
     * The index of the next cluster to be returned
     */
    private int nextCluster;

    public LocsFileReader(final File file) {
        super(file);

        initialize(file);
    }

    public LocsFileReader(final File file, final int lane, final int tile) {
        super(file, lane, tile);

        initialize(file);
    }

    private void initialize(final File file) {
        bbIterator = getFloatIterator(HEADER_SIZE, file);
        final ByteBuffer headerBuf = bbIterator.getHeaderBytes();

        final int firstValue = headerBuf.getInt();
        if (firstValue != BYTES_1_TO_4) {
            throw new IlluminaReaderException("First header byte of locs files should be " + BYTES_1_TO_4 + " value found(" + firstValue + ")");
        }

        final float versionNumber = headerBuf.getFloat();
        if (versionNumber != VERSION) {
            throw new IlluminaReaderException("First header byte of locs files should be " + VERSION + " value found(" + firstValue + ")");
        }

        numClusters = uIntToLong(headerBuf.getInt());
        bbIterator.assertTotalElementsEqual(numClusters * 2);
    }

    @Override
    protected PositionInfo unsafeNextInfo() {
        final float xVal = bbIterator.next();
        final float yVal = bbIterator.next();
        ++nextCluster;
        return new PositionInfo(xVal, yVal, getLane(), getTile());
    }

    @Override
    protected String makeExceptionMsg() {
        return "LocsFileReader(file=" + getFile().getAbsolutePath() + ", numClusters=" + numClusters + ") ";
    }

    @Override
    public boolean hasNext() {
        return nextCluster < numClusters;
    }

    public void close() {
        bbIterator = null;
    }

    public void skipRecords(final int numToSkip) {
        bbIterator.skipElements(numToSkip * 2);
    }
}
