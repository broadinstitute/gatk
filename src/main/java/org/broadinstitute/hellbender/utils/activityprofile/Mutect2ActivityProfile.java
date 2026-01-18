package org.broadinstitute.hellbender.utils.activityprofile;

import htsjdk.samtools.SAMFileHeader;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.engine.spark.AssemblyRegionArgumentCollection;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;
import java.util.stream.IntStream;

public class Mutect2ActivityProfile extends ActivityProfile {

    public Mutect2ActivityProfile(final AssemblyRegionArgumentCollection args, final SAMFileHeader header) {
        // note that in this class maxProbPropagationDistance is interpreted as a maximum phasing distance
        // we don't "propagate" probability like the BandPassActivityFilter
        super(args.maxProbPropagationDistance, args.activeProbThreshold, header);
    }

    @Override
    protected Collection<ActivityProfileState> processState(final ActivityProfileState justAddedState) {
        return Collections.singletonList(justAddedState);
    }

    // if two active loci are separated by this distance or less, they go in the same AssemblyRegion
    // eg if maxPhasingDistance is 10, sites 1 and 11 are assembled together but 1 and 12 are not.
    private int maxPhasingDistance() { return getMaxProbPropagationDistance(); }

    /**
     * Note: when this is called we have already checked readyToPopNextAssemblyRegion().  Therefore, either 1) we are at the
     * end of an interval and have to form a region without knowing about potential variants past the interval's end or 2)
     * the state list spans the maxRegionSize PLUS the maxPhasingDistance.
     *
     * In either case, we have enough information to make regions up to Math.min(stateList.size(), maxRegionSize) bases.
     *
     * Furthermore, any previous region has already checked that it's out of phasing distance with these states.

     */
    @Override
    protected Pair<Integer, Boolean> findSizeOfRegionAndActivity(final int minRegionSize, final int maxRegionSize) {
        // TODO: should this be padded?  The BandPassActivityProfile implicitly pads by propagating
        // TODO: probability with the Gaussian filter
        // TODO: I'm especially nervous about indels -- maybe the actual variant is a base or two before the
        // TODO: active pileup
        Utils.validate(!stateList.isEmpty(), "This code should only be called when there are states with which to form a region.");

        final int maxSize = Math.min(stateList.size(), maxRegionSize);

        // If inactive, extend up to the first active site or the maximum possible size
        // If active, extend until we reach an inactive gap longer than the max phasing distance
        if (getProb(0) < activeProbThreshold) {
            final int numInactiveStates = IntStream.range(0, maxSize)
                    .filter(n -> getProb(n) >= activeProbThreshold)
                    .findFirst().orElse(maxSize);

            return Pair.of(numInactiveStates, false);
        } else {
            int lastActiveIdx = 0;
            int startOfLargestGap = 0;  // the active site right before the longest inactive gap
            int maxGapLength = 0;
            for (int idx = 0; idx < maxSize; idx++) {
                final int gapLength = idx - lastActiveIdx - 1;
                if (gapLength > maxGapLength) {
                    startOfLargestGap = lastActiveIdx;
                    maxGapLength = gapLength;
                }

                if (idx > lastActiveIdx + maxPhasingDistance()) {
                    break;
                } else if (getProb(idx) >= activeProbThreshold) {
                    lastActiveIdx = idx;
                }
            }

            return Pair.of(startOfLargestGap, true);
        }
    }
}
