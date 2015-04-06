package org.broadinstitute.hellbender.tools.picard.illumina.parser;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.StringUtil;


import java.io.File;
import java.util.*;

import static htsjdk.samtools.util.StringUtil.join;

/**
 * Abstract base class for Parsers that open a single tile file at a time and iterate through them.
 */
public abstract class PerTileParser<ILLUMINA_DATA extends IlluminaData> implements IlluminaParser<ILLUMINA_DATA> {
    private final IlluminaFileMap tileToFiles;
    private CloseableIterator<ILLUMINA_DATA> currentIterator;
    private Integer nextTile;
    private Integer currentTile;

    /**
     * Factory method for the iterator of each tile
     */
    protected abstract CloseableIterator<ILLUMINA_DATA> makeTileIterator(final File nextTileFile);

    public PerTileParser(final IlluminaFileMap tilesToFiles) {
        this.tileToFiles = tilesToFiles;
        this.nextTile = tilesToFiles.firstKey();
        this.currentTile = null;
    }

    public PerTileParser(final IlluminaFileMap tilesToFiles, final int nextTile) {
        this.tileToFiles = tilesToFiles;
        this.currentTile = null;
        this.nextTile = nextTile;

        if (!tilesToFiles.containsKey(nextTile)) {
            throw new IllegalArgumentException("NextTile (" + nextTile + ") is not contained by tilesToFiles (" + join(",", new ArrayList<Integer>(tilesToFiles.keySet())));
        }
    }

    /**
     * Return the tile of the NEXT ILLUMINA_DATA object to be returned by the method next.  This might force us to advance to the
     * next file (as it will contains the data for the next) tile/ILLUMINA_DATA object.
     *
     * @return tile number for the next ILLUMINA_DATA object to be returned
     */
    public int getTileOfNextCluster() {
        maybeAdvance();
        return currentTile;
    }

    private void advanceTile() {
        if (nextTile == null) {
            throw new NoSuchElementException("No more tiles to advance!");
        }

        if (currentIterator != null) {
            currentIterator.close();
        }

        currentIterator = makeTileIterator(tileToFiles.get(nextTile));
        currentTile = nextTile;
        nextTile = tileToFiles.higherKey(nextTile);
    }

    public void seekToTile(int oneBasedTileNumber) {
        nextTile = oneBasedTileNumber;

        if (!tileToFiles.containsKey(oneBasedTileNumber)) {
            throw new IlluminaParserException("PerTileParser does not contain key(" + oneBasedTileNumber + ") keys available ("
                    + join(",", new ArrayList<Integer>(tileToFiles.keySet())) + ")");
        }

        if (currentIterator != null) {
            currentIterator.close();
        }
        currentIterator = null;
    }

    public void maybeAdvance() {
        if (!hasNext()) {
            throw new NoSuchElementException();
        }

        if (currentIterator == null || !currentIterator.hasNext()) {
            advanceTile();
        }
    }

    public ILLUMINA_DATA next() {
        maybeAdvance();

        return currentIterator.next();
    }

    public void remove() {
        throw new UnsupportedOperationException();
    }

    public boolean hasNext() {
        // Skip over empty tiles
        while ((currentIterator == null || !currentIterator.hasNext()) && nextTile != null) {
            advanceTile();
        }
        return currentIterator != null && currentIterator.hasNext();
    }

    public void close() {
        if (currentIterator != null) {
            currentIterator.close();
        }
    }

    public void verifyData(List<Integer> tiles, final int[] cycles) {
        final List<Integer> mapTiles = new ArrayList<Integer>(this.tileToFiles.keySet());
        if (!mapTiles.containsAll(tiles)) {
            throw new IlluminaParserException("Missing tiles in PerTileParser expected(" + join(",", tiles)
                    + ") but found (" + join(",", mapTiles) + ")");
        }

        if (!tiles.containsAll(mapTiles)) {
            throw new IlluminaParserException("Extra tiles where found in PerTileParser  expected("
                    + join(",", tiles) + ") but found (" + join(",", mapTiles) + ")");
        }
    }
}
