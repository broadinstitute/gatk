package org.broadinstitute.hellbender.utils.activityprofile;

import htsjdk.samtools.SAMFileHeader;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.engine.spark.AssemblyRegionArgumentCollection;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Collection;
import java.util.Collections;
import java.util.Random;
import java.util.stream.IntStream;

/**
 * Somatic-specific alternative to the default BandPassActivityProfile.  It is based on the following principles:
 *
 *  1) Two events separated by the maximum phasing distance or less should be put in the same active region in
 *      order to be assembled together.
 *  2) Germline variants are not inherently events of interest BUT they should be included in
 *     the active region of nearby somatic variants.
 *
 *  For example, suppose the phasing distance is 100, there are somatic variants at loci 200, 350, and 600; and
 *  there are germline variants at loci 280, 650, and 800.  The somatic variants at 200 and 350 will be in the same
 *  active region along with the germline variant at 280; the somatic variant at 600 will be grouped with the
 *  germline variant at 650; and the germline variant at 800 will not be in an active region.
 */
public class Mutect2ActivityProfile extends ActivityProfile {

    private static final int PADDING = 15;

    private final AssemblyRegionArgumentCollection args;
    private final Random rng = Utils.getRandomGenerator();

    public Mutect2ActivityProfile(final AssemblyRegionArgumentCollection args, final SAMFileHeader header) {
        // In this class maxProbPropagationDistance is interpreted as a maximum phasing distance.
        // It doesn't "propagate" probability like the BandPassActivityFilter
        super(args.maxProbPropagationDistance, args.activeProbThreshold, header);
        this.args = args;
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
    protected Pair<Integer, Boolean> findSizeAndActivityOfRegion(final int minRegionSize, final int maxRegionSize) {
        Utils.validate(!stateList.isEmpty(), "This code should only be called when there are states with which to form a region.");

        final int maxSize = Math.min(stateList.size(), maxRegionSize + PADDING);

        final int numInactiveStates = IntStream.range(0, maxSize)
                .filter(n -> getProb(n) >= activeProbThreshold)
                .findFirst().orElse(maxSize);

        final boolean inactive = numInactiveStates > PADDING;

        // If inactive, extend up to the first active site or the maximum possible size
        // If active, extend until we reach an inactive gap longer than the max phasing distance
        if (inactive) {
            return Pair.of(numInactiveStates - PADDING, false);
        } else {
            boolean isSomatic = false;
            int lastActiveIdx = numInactiveStates;
            int startOfLargestGap = numInactiveStates;  // the active site right before the longest inactive gap
            int maxGapLength = 0;
            for (int idx = numInactiveStates; idx < maxSize; idx++) {
                final int gapLength = idx - lastActiveIdx - 1;
                if (gapLength > maxGapLength) {
                    startOfLargestGap = lastActiveIdx;
                    maxGapLength = gapLength;
                }

                if (idx > lastActiveIdx + maxPhasingDistance()) {
                    break;
                } else if (getProb(idx) >= activeProbThreshold) {
                    isSomatic |= stateList.get(idx).getResultState() == ActivityProfileState.Type.SOMATIC;
                    lastActiveIdx = idx;
                }
            }

            final boolean keepGermline = rng.nextDouble() < args.genotypeGermlineSitesFraction;
            return Pair.of(Math.min(startOfLargestGap + 1 + PADDING, maxSize), isSomatic || keepGermline);
        }
    }
}