package org.broadinstitute.hellbender.tools.spark.pathseq;

import java.util.ArrayList;
import java.util.Collection;

/**
 * Stores taxonomic IDs that were hits of a read pair. ID's can be repeated.
 */
public final class PSHitInfo {
    public final ArrayList<String> taxIDs;
    public final int numMates;

    public PSHitInfo(final Collection<String> taxIDs, final int numMates) {
        this.taxIDs = new ArrayList<>(taxIDs);
        this.numMates = numMates;
    }
}
