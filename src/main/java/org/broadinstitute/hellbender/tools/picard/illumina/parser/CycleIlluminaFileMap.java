package org.broadinstitute.hellbender.tools.picard.illumina.parser;

import java.util.*;

/**
 * For per cycle files.  Maps a Cycle -> Tile -> List<File>
 *
 * @author jburke@broadinstitute.org
 */
class CycleIlluminaFileMap extends TreeMap<Integer, IlluminaFileMap> {
    /**
     * Return a CycleIlluminaFileMap with only the tiles listed and all of the cycles provided.
     * Important NOTE: this DOES NOT eliminate cycles from the cycles parameter passed in that are missing in the cyclesFileIterator of any given lane in the CycleIlluminaFileMap
     */
    public CycleIlluminaFileMap keep(final List<Integer> tilesToKeep, final Set<Integer> cycles) {
        final CycleIlluminaFileMap ciMap = new CycleIlluminaFileMap();

        if (cycles != null) {
            for (final int cycle : cycles) {
                final IlluminaFileMap template = this.get(cycle);
                if (template != null) {
                    ciMap.put(cycle, template.keep(tilesToKeep));
                }
            }
        }

        return ciMap;
    }

    /**
     * Assert that this map has an iterator for all of the expectedTiles and each iterator has expectedCycles number
     * of files.  Also, assert that each cycle file for a given tile is the same size
     *
     * @param expectedTiles  A list of tiles that should be in this map
     * @param expectedCycles The total number of files(cycles) that should be in each CycledFilesIterator
     */
    public void assertValid(final List<Integer> expectedTiles, final int[] expectedCycles) {
        if (size() != expectedCycles.length) {
            throw new IlluminaParserException("Expected CycledIlluminaFileMap to contain " + expectedCycles.length + " cycles but only " + size() + " were found!");
        }
        if (this.firstEntry().getValue().size() != expectedTiles.size()) {
            throw new IlluminaParserException("Expected CycledIlluminaFileMap to contain " + expectedTiles.size()
                    + " tiles but only " + this.firstEntry().getValue().size() + " were found!");
        }
    }

}
