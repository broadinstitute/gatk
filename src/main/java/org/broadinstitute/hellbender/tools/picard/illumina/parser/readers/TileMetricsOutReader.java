package org.broadinstitute.hellbender.tools.picard.illumina.parser.readers;

import org.broadinstitute.hellbender.utils.UnsignedTypeUtil;

import java.io.File;
import java.nio.ByteBuffer;
import java.util.*;

/**
 * Reads a TileMetricsOut file commonly found in the InterOp directory of an Illumina Run Folder.  This
 * reader DOES NOT try to interpret the metrics code or metrics value but instead returns them in what
 * is essentially a struct.
 *
 * File Format:
 * byte 0 (unsigned byte) = The version number which MUST be 2 or an exception will be thrown
 * byte 1 (unsigned byte) = The record size which must be 10 or an exception will be thrown
 * bytes 3 + (current_record * 10) to (current_record * 10 + 10) (TileMetrics Record) = The actual records each of size 10 that
 *          get converted into IlluminaPhasingMetrics objects
 *
 * TileMetrics Record Format:
 * Each 10 byte record is of the following format:
 * byte 0-1 (unsigned short) = lane number
 * byte 2-3 (unsigned short) = tile number
 * byte 4-5 (unisgned short) = metrics code, see Theory of RTA document by Illumina for definition
 * byte 6-9 (float)          = metrics value, see Theory of RTA document by Illumina for definition
 */
public final class TileMetricsOutReader implements Iterator<TileMetricsOutReader.IlluminaTileMetrics> {
    private static final int HEADER_SIZE = 2;
    private static final int EXPECTED_RECORD_SIZE = 10;
    private static final int EXPECTED_VERSION = 2;

    private final BinaryFileIterator<ByteBuffer> bbIterator;

    /**
     * Return a TileMetricsOutReader for the specified file
     * @param tileMetricsOutFile The file to read
     */
    public TileMetricsOutReader(final File tileMetricsOutFile) {
        bbIterator = MMapBackedIteratorFactory.getByteBufferIterator(HEADER_SIZE, EXPECTED_RECORD_SIZE, tileMetricsOutFile);

        final ByteBuffer header = bbIterator.getHeaderBytes();

        //Get the version, should be EXPECTED_VERSION, which is 2
        final int actualVersion = UnsignedTypeUtil.uByteToInt(header.get());
        if(actualVersion != EXPECTED_VERSION) {
            throw new IlluminaReaderException("TileMetricsOutReader expects the version number to be " + EXPECTED_VERSION
                    + ".  Actual Version in Header( " + actualVersion + ")" );
        }

        final int actualRecordSize = UnsignedTypeUtil.uByteToInt(header.get());
        if(EXPECTED_RECORD_SIZE != actualRecordSize) {
            throw new IlluminaReaderException("TileMetricsOutReader expects the record size to be " + EXPECTED_RECORD_SIZE
                    + ".  Actual Record Size in Header( " + actualRecordSize + ")" );
        }
    }

    public boolean hasNext() {
        return bbIterator.hasNext();
    }

    public IlluminaTileMetrics next() {
        if(!hasNext()) {
            throw new NoSuchElementException();
        }
        return new IlluminaTileMetrics(bbIterator.next());
    }

    public void remove() {
        throw new UnsupportedOperationException();
    }

    /**
     * IlluminaPhasingMetrics corresponds to a single record in a TileMetricsOut file
     */
    public static class IlluminaTileMetrics {
        private final IlluminaLaneTileCode laneTileCode;
        private final float metricValue;

        public IlluminaTileMetrics(final ByteBuffer bb) {
            this(UnsignedTypeUtil.uShortToInt(bb.getShort()), UnsignedTypeUtil.uShortToInt(bb.getShort()),
                    UnsignedTypeUtil.uShortToInt(bb.getShort()), bb.getFloat());
        }

        public IlluminaTileMetrics(final int laneNumber, final int tileNumber, final int metricCode, final float metricValue) {
            this.laneTileCode = new IlluminaLaneTileCode(laneNumber, tileNumber, metricCode);
            this.metricValue = metricValue;
        }

        public int getLaneNumber() {
            return laneTileCode.getLaneNumber();
        }

        public int getTileNumber() {
            return laneTileCode.getTileNumber();
        }

        public int getMetricCode() {
            return laneTileCode.getMetricCode();
        }

        public float getMetricValue() {
            return metricValue;
        }

        public IlluminaLaneTileCode getLaneTileCode() {
            return laneTileCode;
        }

        @Override
        public boolean equals(final Object o) {
            if (o instanceof IlluminaTileMetrics) {
                final IlluminaTileMetrics that = (IlluminaTileMetrics) o;
                return laneTileCode == that.laneTileCode && metricValue == that.metricValue; // Identical tile data should render exactly the same float.
            } else {
                return false;
            }
        }
        
        @Override
        public int hashCode() {
            return String.format("%s:%s:%s:%s", laneTileCode.getLaneNumber(), laneTileCode.getTileNumber(), laneTileCode.getMetricCode(), metricValue).hashCode(); // Slow but adequate.
        }
    }

    /** Helper class which captures the combination of a lane, tile & metric code */
    public static class IlluminaLaneTileCode {
        private final int laneNumber;
        private final int tileNumber;
        private final int metricCode;

        public IlluminaLaneTileCode(final int laneNumber, final int tileNumber, final int metricCode) {
            this.laneNumber = laneNumber;
            this.tileNumber = tileNumber;
            this.metricCode = metricCode;
        }

        public int getLaneNumber() {
            return laneNumber;
        }

        public int getTileNumber() {
            return tileNumber;
        }

        public int getMetricCode() {
            return metricCode;
        }

        @Override
        public boolean equals(final Object o) {
            if (o instanceof IlluminaLaneTileCode) {
                final IlluminaLaneTileCode that = (IlluminaLaneTileCode) o;
                return laneNumber == that.laneNumber && tileNumber == that.tileNumber && metricCode == that.metricCode;
            } else {
                return false;
            }
        }

        @Override
        public int hashCode() {
            int result = laneNumber;
            result = 31 * result + tileNumber;
            result = 31 * result + metricCode;
            return result;
        }
    }
}
