package org.broadinstitute.hellbender.tools.picard.illumina.parser.readers;

import htsjdk.samtools.Defaults;
import htsjdk.samtools.util.BlockCompressedInputStream;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.RuntimeIOException;
import org.broadinstitute.hellbender.tools.picard.illumina.parser.BclData;
import org.broadinstitute.hellbender.tools.picard.illumina.parser.TileIndex;
import org.broadinstitute.hellbender.utils.UnsignedTypeUtil;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.*;
import java.util.zip.GZIPInputStream;

/**
 * BCL Files are base call and quality score binary files containing a (base,quality) pair for successive clusters.
 * The file is structured as followed:
 * Bytes 1-4 : unsigned int numClusters
 * Bytes 5-numClusters + 5 : 1 byte base/quality score
 * <p/>
 * The base/quality scores are organized as follows (with one exception, SEE BELOW):
 * The right 2 most bits (these are the LEAST significant bits) indicate the base, where
 * A=00(0x00), C=01(0x01), G=10(0x02), and T=11(0x03)
 * <p/>
 * The remaining bytes compose the quality score which is an unsigned int.
 * <p/>
 * EXCEPTION: If a byte is entirely 0 (e.g. byteRead == 0) then it is a no call, the base
 * becomes '.' and the Quality becomes 2, the default illumina masking value
 * <p/>
 * (E.g. if we get a value in binary of 10001011 it gets transformed as follows:
 * <p/>
 * Value read: 10001011(0x8B)
 * <p/>
 * Quality     Base
 * <p/>
 * 100010      11
 * 00100010    0x03
 * 0x22        T
 * 34          T
 * <p/>
 * So the output base/quality will be a (T/34)
 */
public class BclReader implements CloseableIterator<BclData> {
    private static final byte BASE_MASK = 0x0003;
    private static final int HEADER_SIZE = 4;
    private static final byte[] BASE_LOOKUP = new byte[]{'A', 'C', 'G', 'T'};

    private final InputStream[] streams;
    private final int[] outputLengths;
    int[] numClustersPerCycle;

    private final BclQualityEvaluationStrategy bclQualityEvaluationStrategy;
    private BclData queue = null;

    public BclReader(final List<File> bclsForOneTile, final int[] outputLengths,
                     final BclQualityEvaluationStrategy bclQualityEvaluationStrategy, final boolean seekable) {
        try {
            this.bclQualityEvaluationStrategy = bclQualityEvaluationStrategy;
            this.outputLengths = outputLengths;

            int cycles = 0;
            for (final int outputLength : outputLengths) {
                cycles += outputLength;
            }
            this.streams = new InputStream[cycles];
            this.numClustersPerCycle = new int[cycles];

            final ByteBuffer byteBuffer = ByteBuffer.allocate(HEADER_SIZE);
            byteBuffer.order(ByteOrder.LITTLE_ENDIAN);

            for (int i = 0; i < cycles; ++i) {
                final File bclFile = bclsForOneTile.get(i);
                if (bclFile == null) {
                    close();
                    throw new RuntimeIOException(String.format("Could not find BCL file for cycle %d", i));
                }
                final String filePath = bclFile.getName();
                final boolean isGzip = filePath.endsWith(".gz");
                final boolean isBgzf = filePath.endsWith(".bgzf");
                final InputStream stream = open(bclFile, seekable, isGzip, isBgzf);
                final int read = stream.read(byteBuffer.array());
                if (read != HEADER_SIZE) {
                    close();
                    throw new RuntimeIOException(String.format("BCL %s has invalid header structure.", bclFile.getAbsoluteFile()));
                }
                numClustersPerCycle[i] = byteBuffer.getInt();
                if (!isBgzf && !isGzip) {
                    assertProperFileStructure(bclFile, numClustersPerCycle[i], stream);
                }
                this.streams[i] = stream;
                byteBuffer.clear();
            }
        } catch (final IOException ioe) {
            throw new RuntimeIOException(ioe);
        }
    }

    public static boolean isGzipped(final File file) {
        return file.getAbsolutePath().endsWith(".gz");
    }

    public static boolean isBlockGzipped(final File file) {
        return file.getAbsolutePath().endsWith(".bgzf");
    }

    public static long getNumberOfClusters(final File file) {
        InputStream stream = null;
        try {
            if (isBlockGzipped(file)) stream = new BlockCompressedInputStream(IOUtil.maybeBufferedSeekableStream(file));
            else if (isGzipped(file)) stream = new GZIPInputStream(IOUtil.maybeBufferInputStream(new FileInputStream(file)));
            else stream = IOUtil.maybeBufferInputStream(new FileInputStream(file));

            return getNumberOfClusters(file.getAbsolutePath(), stream);

        } catch (final IOException ioe) {
            throw new IlluminaReaderException("Could not open file " + file.getAbsolutePath() + " to get its cluster count: " + ioe.getMessage(), ioe);
        } finally {
            CloserUtil.close(stream);
        }
    }

    private static long getNumberOfClusters(final String filePath, final InputStream inputStream) {
        final byte[] header = new byte[HEADER_SIZE];

        try {
            final int headerBytesRead = inputStream.read(header);
            if (headerBytesRead != HEADER_SIZE) {
                throw new IlluminaReaderException("Malformed file, expected header of size " + HEADER_SIZE + " but received " + headerBytesRead);
            }
        } catch (final IOException ioe) {
            throw new IlluminaReaderException("Unable to read header for file (" + filePath + ")", ioe);
        }

        final ByteBuffer headerBuf = ByteBuffer.wrap(header);
        headerBuf.order(ByteOrder.LITTLE_ENDIAN);
        return UnsignedTypeUtil.uIntToLong(headerBuf.getInt());
    }


    public BclReader(final File bclFile, final BclQualityEvaluationStrategy bclQualityEvaluationStrategy, final boolean seekable) {
        try {

            this.outputLengths = new int[]{1};
            this.streams = new InputStream[1];
            this.numClustersPerCycle = new int[]{1};
            this.bclQualityEvaluationStrategy = bclQualityEvaluationStrategy;

            final ByteBuffer byteBuffer = ByteBuffer.allocate(HEADER_SIZE);
            final String filePath = bclFile.getName();
            final boolean isGzip = filePath.endsWith(".gz");
            final boolean isBgzf = filePath.endsWith(".bgzf");
            final InputStream stream = open(bclFile, seekable, isGzip, isBgzf);
            final int read = stream.read(byteBuffer.array());

            if (read != HEADER_SIZE) {
                throw new RuntimeIOException(String.format("BCL %s has invalid header structure.", bclFile.getAbsoluteFile()));
            }

            byteBuffer.order(ByteOrder.LITTLE_ENDIAN);
            this.numClustersPerCycle[0] = byteBuffer.getInt();
            if (!isBgzf && !isGzip) {
                assertProperFileStructure(bclFile, this.numClustersPerCycle[0], stream);
            }
            this.streams[0] = stream;
        } catch (final IOException ioe) {
            throw new IlluminaReaderException("IOException opening file " + bclFile.getAbsoluteFile(), ioe);
        }
    }

    void assertProperFileStructure(final File file, final int numClusters, final InputStream stream) {
        final long elementsInFile = file.length() - HEADER_SIZE;
        if (numClusters != elementsInFile) {
            CloserUtil.close(stream);
            throw new IlluminaReaderException("Expected " + numClusters + " in file but found " + elementsInFile);
        }
    }

    InputStream open(final File file, final boolean seekable, final boolean isGzip, final boolean isBgzf) throws IOException {
        final String filePath = file.getAbsolutePath();

        try {
            // Open up a buffered stream to read from the file and optionally wrap it in a gzip stream
            // if necessary
            if (isBgzf) {
                // Only BlockCompressedInputStreams can seek, and only if they are fed a SeekableStream.
                return new BlockCompressedInputStream(IOUtil.maybeBufferedSeekableStream(file));
            } else if (isGzip) {
                if (seekable) {
                    throw new IllegalArgumentException(
                            String.format("Cannot create a seekable reader for gzip bcl: %s.", filePath)
                    );
                }
                return (IOUtil.maybeBufferInputStream(new GZIPInputStream(new FileInputStream(file), Defaults.BUFFER_SIZE / 2),
                        Defaults.BUFFER_SIZE / 2));
            } else {
                if (seekable) {
                    throw new IllegalArgumentException(
                            String.format("Cannot create a seekable reader for provided bcl: %s.", filePath)
                    );
                }
                return IOUtil.maybeBufferInputStream(new FileInputStream(file));
            }
        } catch (final FileNotFoundException fnfe) {
            throw new IlluminaReaderException("File not found: (" + filePath + ")", fnfe);
        } catch (final IOException ioe) {
            throw new IlluminaReaderException("Error reading file: (" + filePath + ")", ioe);
        }
    }

    public void close() {
        for (final InputStream stream : this.streams) {
            CloserUtil.close(stream);
        }
    }

    @Override
    public boolean hasNext() {
        if (queue == null) {
            advance();
        }
        return queue != null;
    }

    public BclData next() {
        if (queue == null) {
            advance();
        }

        final BclData data = queue;
        queue = null;
        return data;
    }

    @Override
    public void remove() {
        throw new UnsupportedOperationException();
    }

    void advance() {
        int totalCycleCount = 0;
        final BclData data = new BclData(outputLengths);
        for (int read = 0; read < outputLengths.length; read++) {
            for (int cycle = 0; cycle < outputLengths[read]; ++cycle) {
                try {
                    final int readByte = this.streams[totalCycleCount].read();
                    if (readByte == -1) {
                        queue = null;
                        return;
                    }

                    if (readByte == 0) {
                        //NO CALL, don't confuse with an A call
                        data.bases[read][cycle] = (byte) '.';
                        data.qualities[read][cycle] = (byte) 2;
                    } else {
                        data.bases[read][cycle] = BASE_LOOKUP[readByte & BASE_MASK];
                        data.qualities[read][cycle] = bclQualityEvaluationStrategy.reviseAndConditionallyLogQuality((byte) (readByte >>> 2));
                    }
                    totalCycleCount++;
                } catch (final IOException ioe) {
                    throw new RuntimeIOException(ioe);
                }

            }
        }
        this.queue = data;
    }

    public static BclReader makeSeekable(final List<File> files, final BclQualityEvaluationStrategy bclQualityEvaluationStrategy, final int[] outputLengths) {
        return new BclReader(files, outputLengths, bclQualityEvaluationStrategy, true);
    }

    public int seek(final List<File> files, final TileIndex tileIndex, final int currentTile) {
        int count = 0;
        int numClustersInTile = 0;
        for (final InputStream inputStream : streams) {
            final TileIndex.TileIndexRecord tileIndexRecord = tileIndex.findTile(currentTile);
            final BclIndexReader bclIndexReader = new BclIndexReader(files.get(count));
            final long virtualFilePointer = bclIndexReader.get(tileIndexRecord.getZeroBasedTileNumber());
            if (!(inputStream instanceof BlockCompressedInputStream)) {
                throw new UnsupportedOperationException("Seeking only allowed on bzgf");
            } else {
                try {
                    if (tileIndex.getNumTiles() != bclIndexReader.getNumTiles()) {
                        throw new IlluminaReaderException(String.format("%s.getNumTiles(%d) != %s.getNumTiles(%d)",
                                tileIndex.getFile().getAbsolutePath(), tileIndex.getNumTiles(), bclIndexReader.getBciFile().getAbsolutePath(), bclIndexReader.getNumTiles()));
                    }
                    ((BlockCompressedInputStream) inputStream).seek(virtualFilePointer);
                    numClustersInTile = tileIndexRecord.getNumClustersInTile();
                } catch (final IOException e) {
                    throw new IlluminaReaderException("Problem seeking to " + virtualFilePointer, e);
                }
            }
            count++;
        }
        return numClustersInTile;
    }
}

