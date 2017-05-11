package org.broadinstitute.hellbender.tools.spark.pathseq;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

/**
 * Stores taxonomic IDs that were hits of a read pair. ID's can be repeated.
 */
public final class PSPathogenAlignmentHit {
    public final List<String> taxIDs;
    public final int numMates;

    public PSPathogenAlignmentHit(final Collection<String> taxIDs, final int numMates) {
        this.taxIDs = new ArrayList<>(taxIDs);
        this.numMates = numMates;
    }
}
