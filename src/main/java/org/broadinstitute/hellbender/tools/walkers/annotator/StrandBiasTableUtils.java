package org.broadinstitute.hellbender.tools.walkers.annotator;

import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.LogManager;

import java.util.ArrayList;
import java.util.List;

/**
 * A class containing convenience methods used in the strand bias annotation calculations
 */
public final class StrandBiasTableUtils {
    private StrandBiasTableUtils() {}

    private final static Logger logger = LogManager.getLogger(StrandBiasTableUtils.class);

    //For now this is only for 2x2 contingency tables
    private static final int ARRAY_DIM = 2;
    private static final int ARRAY_SIZE = ARRAY_DIM * ARRAY_DIM;

    /**
     * Helper function to turn the FisherStrand table into the SB annotation array
     * @param table the table used by the FisherStrand annotation
     * @return the array used by the per-sample Strand Bias annotation
     */
    public static List<Integer> getContingencyArray( final int[][] table ) {
        if(table.length != ARRAY_DIM || table[0].length != ARRAY_DIM) {
            logger.warn("Expecting a " + ARRAY_DIM + "x" + ARRAY_DIM + " strand bias table.");
            return null;
        }

        final List<Integer> list = new ArrayList<>(ARRAY_SIZE);
        list.add(table[0][0]);
        list.add(table[0][1]);
        list.add(table[1][0]);
        list.add(table[1][1]);
        return list;
    }
}