package org.broadinstitute.hellbender.utils.activityprofile;

import htsjdk.samtools.SAMFileHeader;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.engine.spark.AssemblyRegionArgumentCollection;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Random;
import java.util.stream.Collectors;
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

    /**
     * Since this profile doesn't transform states, its state list holds the states Mutect2Engine::isActive returned, so a
     * region can record its active loci, which make up the cluster that formed it, along with their alt allele types.
     */
    @Override
    protected List<ActivityProfileState> activeStatesToRecord(final List<ActivityProfileState> regionStates) {
        return regionStates.stream().filter(state -> state.isActiveProb() >= activeProbThreshold).collect(Collectors.toList());
    }

    // if two active loci are separated by this distance or less, they go in the same AssemblyRegion
    // eg if maxPhasingDistance is 10, sites 1 and 11 are assembled together but 1 and 12 are not.
    private int maxPhasingDistance() { return maxProbPropagationDistance; }

    /**
     * A region may extend PADDING states past maxRegionSize, and to tell whether its last active state ends a cluster of
     * active states we must see maxPhasingDistance states beyond that.  Overriding this method makes the base class wait
     * for that many states before forming a region, except at the end of an interval.
     */
    @Override
    public int getMaxProbPropagationDistance() {
        return maxPhasingDistance() + PADDING;
    }

    /**
     * Note: when this is called we have already checked readyToPopNextAssemblyRegion().  Therefore, either 1) we are at the
     * end of an interval and have to form a region without knowing about potential variants past the interval's end or 2)
     * the state list spans maxRegionSize + PADDING + maxPhasingDistance states (see getMaxProbPropagationDistance).
     *
     * In either case, we can tell whether a cluster of active states that begins within the first maxRegionSize + PADDING
     * states also ends within them: in case 1 because no more states are coming, and in case 2 because we can see
     * maxPhasingDistance states past each of them.
     *
     * The previous region, if any, ended either at the end of a cluster, in which case these states are beyond the phasing
     * distance of its active states, or where we split a cluster too large for one region, in which case that cluster
     * continues here.  Either way, the region we form here starts a new cluster at its first active state.
     */
    @Override
    protected Pair<Integer, Boolean> findSizeAndActivityOfRegion(final int minRegionSize, final int maxRegionSize) {
        Utils.validate(!stateList.isEmpty(), "This code should only be called when there are states with which to form a region.");

        final int maxSize = Math.min(stateList.size(), maxRegionSize + PADDING);

        final int numInactiveStates = IntStream.range(0, maxSize)
                .filter(this::isActiveState)
                .findFirst().orElse(maxSize);

        // A stretch with no active states, or with more than PADDING inactive states before the first active state, is
        // inactive.  In the latter case we leave PADDING inactive states to pad the next active region.
        final boolean inactive = numInactiveStates == maxSize || numInactiveStates > PADDING;

        // If inactive, extend up to the first active site or the maximum possible size
        // If active, extend until we reach an inactive gap longer than the max phasing distance
        if (inactive) {
            return Pair.of(numInactiveStates > PADDING ? numInactiveStates - PADDING : numInactiveStates, false);
        } else {
            // Find the cluster of active states, each within the max phasing distance of the previous one, that begins at
            // the first active state.  We look past maxSize if necessary, stopping at the cluster's first active state past
            // maxSize, if any.  Except at the end of an interval, the state list extends maxPhasingDistance past maxSize (see
            // getMaxProbPropagationDistance), so running out of states means the interval, and hence the cluster, has ended.
            int lastActiveIdx = numInactiveStates;
            for (int idx = numInactiveStates + 1; idx < stateList.size() && lastActiveIdx < maxSize; idx++) {
                if (idx > lastActiveIdx + maxPhasingDistance()) {
                    break;
                } else if (isActiveState(idx)) {
                    lastActiveIdx = idx;
                }
            }

            // If the cluster fits, the region extends PADDING past its last active state.  Otherwise, split the cluster
            // at its longest gap between active states.  The padding is at most the max phasing distance so that it never
            // reaches the next cluster.
            final boolean clusterFits = lastActiveIdx < maxSize;
            final int lastActiveIdxInRegion = clusterFits ? lastActiveIdx : lastActiveStateBeforeLongestGap(numInactiveStates, lastActiveIdx);
            final int regionSize = Math.min(lastActiveIdxInRegion + 1 + Math.min(PADDING, maxPhasingDistance()), maxSize);

            // The region is active if it contains a somatic state, not if the rest of a split cluster does.  Otherwise, if it
            // contains germline states we keep it with probability given by the germline fraction, which Mutect2 sets to zero
            // unless genotyping germline sites.  Regions with only panel of normals (artifact) states are inactive.
            final boolean isSomatic = IntStream.range(0, regionSize).anyMatch(n -> hasType(n, ActivityProfileState.Type.SOMATIC));
            final boolean isGermline = !isSomatic && IntStream.range(0, regionSize).anyMatch(n -> hasType(n, ActivityProfileState.Type.GERMLINE));
            final boolean keepGermline = isGermline && rng.nextDouble() < args.genotypeGermlineSitesFraction;
            return Pair.of(regionSize, isSomatic || keepGermline);
        }
    }

    /**
     * Find where to split a cluster of active states that extends past maxSize: after the active state before maxSize
     * that is followed by the longest gap until the next active state.  Ties favor later states, which yield larger
     * regions, so that, for example, a run of consecutive active states is split at maxSize.
     *
     * @param firstActiveIdx index of the cluster's first active state
     * @param firstActiveIdxPastMaxSize index of the cluster's first active state at or past maxSize
     * @return index of the last active state to include in the region
     */
    private int lastActiveStateBeforeLongestGap(final int firstActiveIdx, final int firstActiveIdxPastMaxSize) {
        int result = firstActiveIdx;
        int longestGap = -1;
        int previousActiveIdx = firstActiveIdx;
        for (int idx = firstActiveIdx + 1; idx <= firstActiveIdxPastMaxSize; idx++) {
            if (isActiveState(idx)) {
                final int gap = idx - previousActiveIdx - 1;
                if (gap >= longestGap) {
                    result = previousActiveIdx;
                    longestGap = gap;
                }
                previousActiveIdx = idx;
            }
        }
        return result;
    }

    private boolean isActiveState(final int index) {
        return getProb(index) >= activeProbThreshold;
    }

    private boolean hasType(final int index, final ActivityProfileState.Type type) {
        return stateList.get(index).getResultState() == type;
    }
}