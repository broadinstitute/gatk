package org.broadinstitute.hellbender.tools.walkers.mutect;

import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.ArrayList;

/**
 * Making this a interface so that in the future, if anyone wants to write a different concensus caller,
 * they may do so. And, using different consensus set, they may get different set of intervals where
 * the top strand and bottom strand disagree, thereby giving us a different truth set for learning PCR errors.
 */
public interface DuplexConsensusCaller {
    void letsDoIt(final ArrayList<GATKRead> duplicateSet, final ReferenceContext referenceContext, final String umi);
}
