package org.broadinstitute.hellbender.tools.picard.illumina.parser;

import java.io.File;
import java.util.*;

/**
 * For "non-cycle" files (files that have multiple cycles per file).  Maps a Tile -> File
 *
 * @author jburke@broadinstitute.org
 */
class IlluminaFileMap extends TreeMap<Integer, File> {

    /**
     * Return a file map that includes only the tiles listed
     */
    public IlluminaFileMap keep(final List<Integer> tilesToKeep) {
        final IlluminaFileMap fileMap = new IlluminaFileMap();
        for (final Integer tile : tilesToKeep) {
            final File file = this.get(tile);
            if (file != null) {
                fileMap.put(tile, file);
            }
        }
        return fileMap;
    }

    /**
     * Return the List of Files in order starting at the given tile and containing all files with tile numbers greater than startingTile that
     * are within this map
     *
     * @param startingTile The first File in the returned list will correspond to this tile
     * @return A List of files for all tiles >= startingTile that are contained in this FileMap
     */
    public List<File> getFilesStartingAt(final int startingTile) {
        return new ArrayList<>(this.tailMap(startingTile).values());
    }
}
