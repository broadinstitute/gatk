package org.broadinstitute.hellbender.tools.picard.illumina.parser;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;


import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.*;

import static htsjdk.samtools.util.CloserUtil.close;
import static htsjdk.samtools.util.IOUtil.maybeBufferInputStream;
import static java.lang.String.format;
import static java.nio.ByteBuffer.allocate;
import static java.nio.ByteOrder.LITTLE_ENDIAN;

/**
 * Load a file containing 8-byte records like this:
 * tile number: 4-byte int
 * number of clusters in tile: 4-byte int
 * Number of records to read is determined by reaching EOF.
 */
public class TileIndex implements Iterable<TileIndex.TileIndexRecord> {
    private final File tileIndexFile;
    private final List<TileIndexRecord> tiles = new ArrayList<TileIndexRecord>();

    TileIndex(final File tileIndexFile) {
        try {
            this.tileIndexFile = tileIndexFile;
            final InputStream is = maybeBufferInputStream(new FileInputStream(tileIndexFile));
            final ByteBuffer buf = allocate(8);
            buf.order(LITTLE_ENDIAN);
            int absoluteRecordIndex = 0;
            int numTiles = 0;
            while (readTileIndexRecord(buf.array(), buf.capacity(), is)) {
                buf.rewind();
                buf.limit(buf.capacity());
                final int tile = buf.getInt();
                // Note: not handling unsigned ints > 2^31, but could if one of these exceptions is thrown.
                if (tile < 0)
                    throw new IlluminaParserException("Tile number too large in " + tileIndexFile.getAbsolutePath());
                final int numClusters = buf.getInt();
                if (numClusters < 0)
                    throw new IlluminaParserException("Cluster size too large in " + tileIndexFile.getAbsolutePath());
                tiles.add(new TileIndexRecord(tile, numClusters, absoluteRecordIndex, numTiles++));
                absoluteRecordIndex += numClusters;
            }
            close(is);
        } catch (final IOException e) {
            throw new IlluminaParserException("Problem reading " + tileIndexFile.getAbsolutePath(), e);
        }
    }

    public File getFile() {
        return tileIndexFile;
    }

    public int getNumTiles() {
        return tiles.size();
    }

    private boolean readTileIndexRecord(final byte[] buf, final int numBytes, final InputStream is) throws IOException {
        int totalBytesRead = 0;
        while (totalBytesRead < numBytes) {
            final int bytesRead = is.read(buf, totalBytesRead, numBytes - totalBytesRead);
            if (bytesRead == -1) {
                if (totalBytesRead != 0) {
                    throw new IlluminaParserException(tileIndexFile.getAbsolutePath() + " has incomplete last block");
                } else return false;
            }
            totalBytesRead += bytesRead;
        }
        return true;
    }

    public List<Integer> getTiles() {
        final List<Integer> ret = new ArrayList<Integer>(tiles.size());
        for (final TileIndexRecord rec : tiles) ret.add(rec.tile);
        return ret;
    }

    public List<String> verify(final List<Integer> expectedTiles) {
        final Set<Integer> tileSet = new HashSet<Integer>(tiles.size());
        for (final TileIndexRecord rec : tiles) tileSet.add(rec.tile);
        final List<String> failures = new LinkedList<String>();
        for (final int expectedTile : expectedTiles) {
            if (!tileSet.contains(expectedTile)) {
                failures.add("Tile " + expectedTile + " not found in " + tileIndexFile.getAbsolutePath());
            }
        }
        return failures;
    }

    @Override
    public Iterator<TileIndexRecord> iterator() {
        return tiles.iterator();
    }

    /**
     * @throws java.util.NoSuchElementException if tile is not found
     */
    public TileIndexRecord findTile(final int tileNumber) {
        for (final TileIndexRecord rec : this) {
            if (rec.tile == tileNumber) return rec;
            if (rec.tile > tileNumber) {
                break;
            }
        }
        throw new NoSuchElementException(format("Tile %d not found in %s", tileNumber, tileIndexFile));
    }

    public static class TileIndexRecord {
        /**
         * Number of the tile, e.g. 11101.  These don't necessarily start at 0, and there may be gaps.
         */
        final int tile;

        final int numClustersInTile;

        public int getNumClustersInTile() {
            return numClustersInTile;
        }

        public int getZeroBasedTileNumber() {
            return zeroBasedTileNumber;
        }

        /**
         * I.e. the sum of numClustersInTile for all tiles preceding this one.
         */
        final int indexOfFirstClusterInTile;

        /**
         * A contiguous numbering of tiles starting at 0.
         */
        final int zeroBasedTileNumber;

        private TileIndexRecord(final int tile, final int numClustersInTile, final int indexOfFirstClusterInTile, final int zeroBasedTileNumber) {
            this.tile = tile;
            this.numClustersInTile = numClustersInTile;
            this.indexOfFirstClusterInTile = indexOfFirstClusterInTile;
            this.zeroBasedTileNumber = zeroBasedTileNumber;
        }
    }
}
