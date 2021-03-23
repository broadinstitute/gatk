package org.broadinstitute.hellbender.utils.pileup;

import com.google.common.primitives.Bytes;
import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.engine.AlignmentContext;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.LocusWalkerByInterval;
import org.broadinstitute.hellbender.engine.ReferenceContext;

import java.util.List;
import java.util.Set;

public class PileupToHaplotype extends LocusWalkerByInterval {
    @Override
    public List<Locatable> getIntervalObjectsToQueryOver() {
        return null;
    }

    @Override
    public void apply(AlignmentContext alignmentContext, ReferenceContext referenceContext, FeatureContext featureContext, Set<Locatable> activeIntervals) {
        byte refBase = referenceContext.getBase();
        for (Locatable loc : activeIntervals) {
            if (alignmentContext.contains(loc)) {
                List<Byte> bases = Bytes.asList(alignmentContext.getBasePileup().getBases());
                // this contains each base at this site
            }
        }
    }

    @Override
    public void onIntervalStart(Locatable activeInterval) {

    }

    @Override
    public void onIntervalEnd(Locatable activeInterval) {

    }
}
