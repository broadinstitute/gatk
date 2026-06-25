package org.broadinstitute.hellbender.utils.activityprofile;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.engine.AssemblyRegion;
import org.broadinstitute.hellbender.engine.spark.AssemblyRegionArgumentCollection;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;
import java.util.stream.IntStream;

/**
 * Class holding information about per-base activity scores for
 * assembly region traversal
 */
public class ActivityProfile {

    public enum ProfileType {
        BASE, BAND_PASS, MUTECT2
    }

    public static ActivityProfile create(final AssemblyRegionArgumentCollection args, final SAMFileHeader header) {
        return switch (args.activityProfileType) {
            case BASE -> new ActivityProfile(args.maxProbPropagationDistance, args.activeProbThreshold, header);
            case BAND_PASS -> new BandPassActivityProfile(args.maxProbPropagationDistance, args.activeProbThreshold,
                    BandPassActivityProfile.MAX_FILTER_SIZE, BandPassActivityProfile.DEFAULT_SIGMA, header);
            case MUTECT2 -> new Mutect2ActivityProfile(args, header);
        };
    }

    protected final List<ActivityProfileState> stateList;

    protected final int maxProbPropagationDistance;
    protected final double activeProbThreshold;

    protected SimpleInterval regionStartLoc = null;
    protected SimpleInterval regionStopLoc = null;

    protected SAMFileHeader samHeader;

    /**
     * A cached value of the regionStartLoc contig length, to make calls to
     * getCurrentContigLength efficient
     */
    protected int contigLength = -1;

    /**
     * Create a empty ActivityProfile, restricting output to profiles overlapping intervals, if not null
     * @param maxProbPropagationDistance region probability propagation distance beyond its maximum size
     * @param activeProbThreshold threshold for the probability of a profile state being active
     */
    public ActivityProfile(final int maxProbPropagationDistance, final double activeProbThreshold, final SAMFileHeader header) {
        this.stateList = new ArrayList<>();
        this.maxProbPropagationDistance = maxProbPropagationDistance;
        this.activeProbThreshold = activeProbThreshold;
        this.samHeader = header;
    }

    @Override
    public String toString() {
        return "ActivityProfile{" +
                "start=" + regionStartLoc +
                ", stop=" + regionStopLoc +
                '}';
    }

    /**
     * How far away can probability mass be moved around in this profile?
     *
     * This distance puts an upper limit on how far, in bp, we will ever propagate probability mass around
     * when adding a new ActivityProfileState.  For example, if the value of this function is
     * 10, and you are looking at a state at bp 5, and we know that no states beyond 5 + 10 will have
     * their probability propagated back to that state.
     *
     * @return a positive integer distance in bp
     */
    public int getMaxProbPropagationDistance() {
        return maxProbPropagationDistance;
    }

    /**
     * How many profile results are in this profile?
     * @return the number of profile results
     */
    public int size() {
        return stateList.size();
    }

    /**
     * Is this profile empty? (ie., does it contain no ActivityProfileStates?)
     * @return true if the profile is empty (ie., contains no ActivityProfileStates)
     */
    public boolean isEmpty() {
        return stateList.isEmpty();
    }

    /**
     * Get the span of this activity profile, which is from the start of the first state to the stop of the last
     * @return a potentially null SimpleInterval.  Will be null if this profile is empty
     */
    public SimpleInterval getSpan() {
        return isEmpty() ? null : regionStartLoc.spanWith(regionStopLoc);
    }

    public String getContig() {
        return regionStartLoc.getContig();
    }

    public int getEnd() {
        return regionStopLoc.getEnd();
    }

    /**
     * Get the list of activity profile results in this object
     * @return a non-null, ordered list of activity profile results
     */
    protected List<ActivityProfileState> getStateList() {
        return stateList;
    }

    /**
     * Get the probabilities of the states as a single linear array of doubles
     * @return a non-null array
     */
    protected double[] getProbabilitiesAsArray() {
        final double[] probs = new double[getStateList().size()];
        int i = 0;
        for ( final ActivityProfileState state : getStateList() ) {
            probs[i++] = state.isActiveProb();
        }
        return probs;
    }

    /**
     * Helper function that gets the interval for a site offset from relativeLoc, protecting ourselves from
     * falling off the edge of the contig.
     *
     * @param relativeLoc the location offset is relative to
     * @param offset the offset from relativeLoc where we'd like to create a GenomeLoc
     * @return a SimpleInterval with relativeLoc.start + offset, if this is on the contig, null otherwise
     */
    protected SimpleInterval getLocForOffset(final SimpleInterval relativeLoc, final int offset) {
        Utils.nonNull(relativeLoc);

        final int start = relativeLoc.getStart() + offset;
        if ( start < 1 || start > getCurrentContigLength() ) {
            return null;
        } else {
            return new SimpleInterval(regionStartLoc.getContig(), start, start);
        }
    }

    /**
     * Get the length of the current contig
     * @return the length in bp
     */
    private int getCurrentContigLength() {
        return contigLength;
    }

    // --------------------------------------------------------------------------------
    //
    // routines to add states to a profile
    //
    // --------------------------------------------------------------------------------

    /**
     * Add the next ActivityProfileState to this profile.
     *
     * Must be contiguous with the previously added result, or an IllegalArgumentException will be thrown
     *
     * @param state a well-formed ActivityProfileState result to incorporate into this profile
     */
    public void add(final ActivityProfileState state) {
        Utils.nonNull(state);
        final SimpleInterval loc = state.getLoc();

        if ( regionStartLoc == null ) {
            regionStartLoc = loc;
            regionStopLoc = loc;
            contigLength = samHeader.getSequence(regionStartLoc.getContig()).getSequenceLength();
        } else {
            Utils.validateArg( regionStopLoc.getStart() == loc.getStart() - 1, () ->
                    "Bad add call to ActivityProfile: loc " + loc + " not immediately after last loc " + regionStopLoc);
            regionStopLoc = loc;
        }

        final Collection<ActivityProfileState> processedStates = processState(state);
        for ( final ActivityProfileState processedState : processedStates ) {
            incorporateSingleState(processedState);
        }
    }

    /**
     * Incorporate a single activity profile state into the current list of states
     *
     * If state's position occurs immediately after the last position in this profile, then
     * the state is appended to the state list.  If it's within the existing states list,
     * the prob of stateToAdd is added to its corresponding state in the list.  If the
     * position would be before the start of this profile, stateToAdd is simply ignored.
     *
     * @param stateToAdd the state we want to add to the states list
     */
    private void incorporateSingleState(final ActivityProfileState stateToAdd) {
        Utils.nonNull(stateToAdd);
        final int position = stateToAdd.getOffset(regionStartLoc);
        // should we allow this?  probably not
        Utils.validateArg(position <= size(), () -> "Must add state contiguous to existing states: adding " + stateToAdd);

        if ( position >= 0 ) {
            // ignore states starting before this region's start
            if ( position < size() ) {
                stateList.get(position).setIsActiveProb(stateList.get(position).isActiveProb() + stateToAdd.isActiveProb());
            } else {
                Utils.validateArg(position == size(), "position == size but it wasn't");
                stateList.add(stateToAdd);
            }
        }
    }

    /**
     * Process justAddedState, returning a collection of derived states that actually be added to the stateList
     *
     * The purpose of this function is to transform justAddedStates, if needed, into a series of atomic states
     * that we actually want to track.  For example, if state is for soft clips, we transform that single
     * state into a list of states that surround the state up to the distance of the soft clip.
     *
     * Can be overridden by subclasses to transform states in any way
     *
     * There's no particular contract for the output states, except that they can never refer to states
     * beyond the current end of the stateList unless the explicitly include preceding states before
     * the reference.  So for example if the current state list is [1, 2, 3] this function could return
     * [1,2,3,4,5] but not [1,2,3,5].
     *
     * @param justAddedState the state our client provided to use to add to the list
     * @return a list of derived states that should actually be added to this profile's state list
     */
    protected Collection<ActivityProfileState> processState(final ActivityProfileState justAddedState) {
        if ( justAddedState.getResultState().equals(ActivityProfileState.Type.HIGH_QUALITY_SOFT_CLIPS) ) {
            // special code to deal with the problem that high quality soft clipped bases aren't added to pileups
            final List<ActivityProfileState> states = new ArrayList<>();
            // add no more than the max prob propagation distance num HQ clips
            final int numHQClips = Math.min(justAddedState.getResultValue().intValue(), getMaxProbPropagationDistance());
            for( int i = - numHQClips; i <= numHQClips; i++ ) {
                final SimpleInterval loc = getLocForOffset(justAddedState.getLoc(), i);
                if ( loc != null ) {
                    states.add(new ActivityProfileState(loc, justAddedState.isActiveProb()));
                }
            }

            return states;
        } else {
            return Collections.singletonList(justAddedState);
        }
    }

    // --------------------------------------------------------------------------------
    //
    // routines to get active regions from the profile
    //
    // --------------------------------------------------------------------------------

    /**
     * Get the next completed assembly regions from this profile, and remove all states supporting them from this profile
     *
     * Takes the current profile and finds all of the active / inactive from the start of the profile that are
     * ready.  By ready we mean unable to have their probability modified any longer by future additions to the
     * profile.  The regions that are popped off the profile take their states with them, so the start of this
     * profile will always be after the end of the last region returned here.
     *
     * The regions are returned sorted by genomic position.
     *
     * This function may not return anything in the list, if no regions are ready
     *
     * No returned region will be larger than maxRegionSize.
     *
     * @param assemblyRegionExtension the extension value to provide to the constructed regions
     * @param minRegionSize the minimum region size, in the case where we have to cut up regions that are too large
     * @param maxRegionSize the maximize size of the returned region
     * @param forceConversion if true, we'll return a region whose end isn't sufficiently far from the end of the
     *                        stateList.  Used to close out the active region when we've hit some kind of end (such
     *                        as the end of the contig)
     * @return a non-null list of active regions
     */
    public List<AssemblyRegion> popReadyAssemblyRegions( final int assemblyRegionExtension, final int minRegionSize, final int maxRegionSize, final boolean forceConversion ) {
        Utils.validateArg(assemblyRegionExtension >= 0, () -> "assemblyRegionExtension must be >= 0 but got " + assemblyRegionExtension);
        Utils.validateArg( minRegionSize > 0, () -> "minRegionSize must be >= 1 but got " + minRegionSize);
        Utils.validateArg( maxRegionSize > 0, () -> "maxRegionSize must be >= 1 but got " + maxRegionSize);

        final List<AssemblyRegion> regions = new ArrayList<>();

        while (readyToPopNextAssemblyRegion(maxRegionSize, forceConversion)) {
            final AssemblyRegion nextRegion = popNextReadyAssemblyRegion(assemblyRegionExtension, minRegionSize, maxRegionSize, forceConversion);
            regions.add(nextRegion);
        }

        return regions;
    }

    /**
     * Is there a large enough span of states in the state list to find the next assembly region, or are we
     * at the end of an interval and forcing conversion?
     */
    private boolean readyToPopNextAssemblyRegion(final int maxRegionSize, final boolean forceConversion) {
        return !stateList.isEmpty() && (forceConversion || stateList.size() >= maxRegionSize + getMaxProbPropagationDistance());
    }

    /**
     * Helper function for popReadyActiveRegions that pops the first ready region off the front of this profile
     *
     * If a region is returned, modifies the state of this profile so that states used to make the region are
     * no longer part of the profile.  Associated information (like the region start position) of this profile
     * are also updated.
     *
     * @param assemblyRegionExtension the extension value to provide to the constructed regions
     * @param minRegionSize the minimum region size, in the case where we have to cut up regions that are too large
     * @param maxRegionSize the maximize size of the returned region
     * @param forceConversion if true, we'll return a region whose end isn't sufficiently far from the end of the
     *                        stateList.  Used to close out the active region when we've hit some kind of end (such
     *                        as the end of the contig)
     * @return a fully formed assembly region, or null if none can be made
     */
    private AssemblyRegion popNextReadyAssemblyRegion( final int assemblyRegionExtension, final int minRegionSize, final int maxRegionSize, final boolean forceConversion ) {
        // If we are flushing the activity profile we need to trim off the excess states so that we don't create regions outside of our current processing interval
        if( forceConversion ) {
            stateList.subList(getSpan().size(), stateList.size()).clear();
        }

        final ActivityProfileState first = stateList.get(0);
        final boolean isActiveRegion = first.isActiveProb() > activeProbThreshold;
        final int sizeOfNextRegion = findSizeOfRegion(minRegionSize, maxRegionSize);

        // we need to create the active region, and clip out the states we're extracting from this profile
        final List<ActivityProfileState> sub = stateList.subList(0, sizeOfNextRegion);
        sub.clear();

        // update the start and stop locations as necessary
        if ( stateList.isEmpty() ) {
            regionStartLoc = regionStopLoc = null;
        } else {
            regionStartLoc = stateList.get(0).getLoc();
        }
        final SimpleInterval regionLoc = new SimpleInterval(first.getLoc().getContig(), first.getLoc().getStart(), first.getLoc().getStart() + sizeOfNextRegion - 1);
        return new AssemblyRegion(regionLoc, isActiveRegion, assemblyRegionExtension, samHeader);
    }

    /**
     * Find the size of the current region.  If the first element in the stateList is below the probability threshold,
     * return the number of consecutive states below the threshold.  Otherwise return the number of consecutive states
     * above the probability threshold, unless this yields too large an active region, in which case we split the active
     * region at a local probability minimum.
     *
     * @param minRegionSize the minimum region size, in the case where we have to cut up regions that are too large
     * @param maxRegionSize the maximize size of the returned region
     * @return the size of this region
     */
    private int findSizeOfRegion(final int minRegionSize, final int maxRegionSize) {
        Utils.validate(!stateList.isEmpty(), "state list is empty");
        Utils.validateArg(minRegionSize >= 1, "minRegionSize must be >= 1");
        final boolean aboveThreshold = getProb(0) > activeProbThreshold;
        final int maxContiguousSize = Math.min(stateList.size(), maxRegionSize);
        final int contiguousSize = IntStream.range(0, maxContiguousSize)
                .filter(n -> getProb(n) > activeProbThreshold != aboveThreshold)
                .findFirst().orElse(maxContiguousSize);

        final boolean needToSplit = aboveThreshold && contiguousSize == maxRegionSize;
        // If the active region is too big, split at the lowest local minimum, ties favoring later cuts, defaulting to
        // the max region size if no local minima exist
        if (needToSplit) {
            return IntStream.range(minRegionSize - 1, contiguousSize).filter(this::isLocalMinimum).boxed()
                    .sorted(Comparator.comparingDouble((Integer n) -> getProb(n)).thenComparingInt(n -> -n))
                    .findFirst().orElse(contiguousSize - 1) + 1;
        } else {
            return contiguousSize;
        }
    }

    /**
     * Helper function to get the probability of the state at offset index
     * @param index a valid offset into the state list
     * @return the isActiveProb of the state at index
     */
    private double getProb(final int index) {
        Utils.validIndex(index, stateList.size());

        return stateList.get(index).isActiveProb();
    }

    /**
     * Is the probability at index in a local minimum?
     *
     * Checks that the probability at index is <= both the probabilities to either side.
     * Returns false if index is at the end or the start of the state list.
     *
     * @param index the index of the state we want to test
     * @return true if prob at state is a minimum, false otherwise
     */
    private boolean isLocalMinimum(final int index) {
        return index > 0 && index < stateList.size() - 1 && getProb(index) <= getProb(index + 1) &&
                getProb(index) < getProb(index - 1);
    }
}
