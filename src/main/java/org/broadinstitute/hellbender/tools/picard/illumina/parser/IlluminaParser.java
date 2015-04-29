package org.broadinstitute.hellbender.tools.picard.illumina.parser;

import java.util.*;

/**
 * Interface for classes that parse information out of the Illumina Pipeline
 *
 * @author jburke@broadinstitute.org
 */
interface IlluminaParser<DATA_TYPE extends IlluminaData> extends Iterator<DATA_TYPE> {
    /**
     * Jump so that the next record returned will be from the specified tile.
     */
    void seekToTile(int oneBasedTileNumber);

    /**
     * Read the next read's set of data and set it into the provided data object.  The object must have
     * the appropriate IlluminaEndData objects set into it for first end, second end, barcode.
     */
    DATA_TYPE next();

    /**
     * Is there a DATA_TYPE object for another cluster remaining.
     *
     * @return TRUE if there is a DATA_TYPE object for the next cluster that can be provided by
     * next
     */
    boolean hasNext();

    /**
     * Get the tile for the NEXT DATA_TYPE object that will be returned by this parser.  This should
     * be called BEFORE next if you want the tile for the value returned by next
     */
    public int getTileOfNextCluster();

    /**
     * Given the expected tiles and cycles for this run, make sure this parser can provide data for
     * all tiles/cycles or throws an exception if it's missing any required data or data structures
     * it relies on do not disagree with the provided tiles/cycles
     *
     * @param tiles  The number of tiles in the current run
     * @param cycles The number of cycles in the current run
     */
    void verifyData(final List<Integer> tiles, final int[] cycles);

    void close();

}
