package org.broadinstitute.hellbender.tools.picard.illumina.parser;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.StringUtil;

import org.broadinstitute.hellbender.utils.UnsignedTypeUtil;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.MappedByteBuffer;
import java.nio.channels.FileChannel;
import java.util.*;

import static htsjdk.samtools.util.CloserUtil.close;
import static htsjdk.samtools.util.StringUtil.bytesToString;
import static htsjdk.samtools.util.StringUtil.stringToBytes;
import static java.nio.ByteBuffer.allocate;
import static java.nio.ByteOrder.LITTLE_ENDIAN;
import static java.nio.channels.FileChannel.MapMode;
import static java.nio.channels.FileChannel.MapMode.READ_ONLY;
import static org.broadinstitute.hellbender.tools.picard.illumina.parser.IntensityChannel.values;
import static org.broadinstitute.hellbender.utils.UnsignedTypeUtil.uShortToInt;

/**
 * Read a .cnf (binary noise) or .cif (binary intensity) file.  A file in this format contains
 * 1 or more cycles of data for a set of clusters, with 4 values per cycle, one for each channel.
 * A file can store its values in either a byte or a short per value, but the API treats them all as shorts.
 * This class does not distinguish btw CIF and CNF files.
 *
 * @author jburke@broadinstitute.org
 */
final class ClusterIntensityFileReader {

    private static final byte[] IDENTIFIER = stringToBytes("CIF");
    private static final byte FILE_VERSION = 1;
    private static final int HEADER_SIZE = 13;
    private static final int NUM_CHANNELS = values().length;

    // Just for error reporting
    private final File file;

    /**
     * The entire file is mmapped
     */
    private final MappedByteBuffer buf;
    private final ClusterIntensityFileHeader header;

    // Precomputed for speed, I hope.
    private final int cycleSize;
    private final int channelSize;

    public static class ClusterIntensityFileHeader {
        public final int elementSize;
        public final int firstCycle;
        public final int numCycles;
        public final int numClusters;

        public ClusterIntensityFileHeader(final byte[] headerBytes, final File file) {
            if (headerBytes.length < HEADER_SIZE) {
                throw new IlluminaParserException("Bytes past to header constructor are too short excpected(" + HEADER_SIZE + ") received (" + headerBytes.length);
            }

            ByteBuffer buf = allocate(headerBytes.length); //for doing some byte conversions
            buf.order(LITTLE_ENDIAN);
            buf.put(headerBytes);
            buf.position(0);

            final byte[] identifierBuf = new byte[IDENTIFIER.length];
            buf.get(identifierBuf);
            if (!Arrays.equals(identifierBuf, IDENTIFIER)) {
                throw new IlluminaParserException("Cluster intensity file " + file + " contains unexpected header: " +
                        bytesToString(identifierBuf));
            }
            final byte fileVersion = buf.get();
            if (fileVersion != FILE_VERSION) {
                throw new IlluminaParserException("Cluster intensity file " + file + " contains unexpected version: " + fileVersion);
            }
            elementSize = buf.get();
            if (elementSize < 1 || elementSize > 2) {
                throw new IlluminaParserException("Cluster intensity file " + file + " contains unexpected element size: " + elementSize);
            }
            // convert these to unsigned
            firstCycle = uShortToInt(buf.getShort());
            numCycles = uShortToInt(buf.getShort());
            if (numCycles == 0) {
                throw new IlluminaParserException("Cluster intensity file " + file + " has zero cycles.");
            }
            numClusters = buf.getInt();
            if (numClusters < 0) {
                // It is possible for there to be no clusters in a tile.
                throw new IlluminaParserException("Cluster intensity file " + file + " has negative number of clusters: " + numClusters);
            }
        }
    }

    /**
     * Prepare to parse a CIF or CNF file.
     *
     * @param file The file to be parsed.
     */
    public ClusterIntensityFileReader(final File file) {
        try {
            this.file = file;
            final FileInputStream is = new FileInputStream(this.file);
            final FileChannel channel = is.getChannel();
            final long fileSize = channel.size();
            buf = channel.map(READ_ONLY, 0, fileSize);
            buf.order(LITTLE_ENDIAN);
            close(channel);
            close(is);
            final byte[] headerBytes = new byte[HEADER_SIZE];
            buf.get(headerBytes);
            this.header = new ClusterIntensityFileHeader(headerBytes, this.file);
        } catch (IOException e) {
            throw new RuntimeIOException("IOException opening cluster intensity file " + file, e);
        }
        cycleSize = NUM_CHANNELS * header.numClusters * header.elementSize;
        channelSize = header.numClusters * header.elementSize;
    }

    /**
     * Get the value for the given args.  Value is returned as a signed short regardless of whether storage is
     * in bytes or shorts.
     *
     * @param cluster 0-based cluster number.
     * @param channel Which channel is desired.
     * @param cycle   Absolute cycle number.  E.g. if the first cycle in the file is N, then the first value that can
     *                be fetched is cycle=N
     * @return Intensity or noise (depending on whether this is a CIF or CNF file).
     */
    public short getValue(final int cluster, final IntensityChannel channel, final int cycle) {
        if (cycle < header.firstCycle || cycle >= header.firstCycle + header.numCycles) {
            throw new IllegalArgumentException("Requested cycle (" + cycle + ") number out of range.  First cycle=" +
                    header.firstCycle + "; numCycles=" + header.numCycles);
        }
        if (cluster < 0 || cluster >= header.numClusters) {
            throw new IllegalArgumentException("Requested cluster (" + cluster + ") number out of range. numClustersInTile=" + header.numClusters);
        }
        final int relativeCycle = cycle - header.firstCycle;
        final int position = HEADER_SIZE + relativeCycle * cycleSize + channel.ordinal() * channelSize + cluster * header.elementSize;
        buf.position(position);
        if (header.elementSize == 1) {
            return buf.get();
        } else {
            return buf.getShort();
        }
    }

    public File getFile() {
        return file;
    }

    /**
     * @return The first (one-based) cycle stored in this file.
     */
    public int getFirstCycle() {
        return header.firstCycle;
    }

    /**
     * @return Number of clusters stored in this file.
     */
    public int getNumClusters() {
        return header.numClusters;
    }

    /**
     * @return Number of cycles stored in this file.
     */
    public int getNumCycles() {
        return header.numCycles;
    }

    /**
     * @return the size of one intensity value for one channel in this file.
     */
    public int getElementSize() {
        return header.elementSize;
    }

    public static ClusterIntensityFileHeader readHeaders(final File intensityFile) {
        FileInputStream reader = null;
        byte[] headerBytes = new byte[HEADER_SIZE];
        int bytesRead = 0;
        try {
            reader = new FileInputStream(intensityFile);
            bytesRead = reader.read(headerBytes);
        } catch (FileNotFoundException fnfExc) {
            throw new IlluminaParserException("Error opening intensity file (" + intensityFile.getAbsolutePath() + ")", fnfExc);
        } catch (IOException ioExc) {
            throw new IlluminaParserException("Error reading values from header for intensity file (" + intensityFile.getAbsolutePath() + ")", ioExc);
        } finally {
            close(reader);
        }

        if (bytesRead != HEADER_SIZE)
            throw new IlluminaParserException("Error reading intensity file header, too few bytes read, expected( " + HEADER_SIZE + ") read(" + bytesRead + ")");

        return new ClusterIntensityFileHeader(headerBytes, intensityFile);
    }
}
