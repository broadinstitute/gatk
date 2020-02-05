package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.SAMFileHeader;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.spark.AssemblyRegionWalkerSpark;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.activityprofile.ActivityProfile;
import org.broadinstitute.hellbender.utils.activityprofile.ActivityProfileState;
import org.broadinstitute.hellbender.utils.activityprofile.BandPassActivityProfile;

import java.util.*;

/**
 * Given an iterator of {@link ActivityProfileState}, finds {@link AssemblyRegion}s.
 *
 * This iterator and {@link ActivityProfileStateIterator} represent the core of the {@link AssemblyRegionWalkerSpark} traversal.
 */
public class AssemblyRegionFromActivityProfileStateIterator implements Iterator<AssemblyRegion> {
    private static final Logger logger = LogManager.getLogger(AssemblyRegionFromActivityProfileStateIterator.class);

    private final SAMFileHeader readHeader;
    private final int minRegionSize;
    private final int maxRegionSize;
    private final int assemblyRegionPadding;

    private AssemblyRegion readyRegion;
    private Queue<AssemblyRegion> pendingRegions;
    private final Iterator<ActivityProfileState> activityProfileStateIterator;
    private final ActivityProfile activityProfile;

    /**
     * Constructs an AssemblyRegionIterator over a provided read shard
     *
     * @param readHeader header for the reads
     * @param minRegionSize minimum size of an assembly region
     * @param maxRegionSize maximum size of an assembly region
     * @param assemblyRegionPadding number of bases of padding on either side of an assembly region
     * @param activeProbThreshold minimum probability for a locus to be considered active
     * @param maxProbPropagationDistance upper limit on how many bases away probability mass can be moved around
     *                                   when calculating the boundaries between active and inactive assembly regions
     */
    public AssemblyRegionFromActivityProfileStateIterator(final Iterator<ActivityProfileState> activityProfileStateIterator,
                                                          final SAMFileHeader readHeader,
                                                          final int minRegionSize,
                                                          final int maxRegionSize,
                                                          final int assemblyRegionPadding,
                                                          final double activeProbThreshold,
                                                          final int maxProbPropagationDistance) {

        Utils.nonNull(readHeader);
        Utils.validateArg(minRegionSize >= 1, "minRegionSize must be >= 1");
        Utils.validateArg(maxRegionSize >= 1, "maxRegionSize must be >= 1");
        Utils.validateArg(minRegionSize <= maxRegionSize, "minRegionSize must be <= maxRegionSize");
        Utils.validateArg(assemblyRegionPadding >= 0, "assemblyRegionPadding must be >= 0");
        Utils.validateArg(activeProbThreshold >= 0.0, "activeProbThreshold must be >= 0.0");
        Utils.validateArg(maxProbPropagationDistance >= 0, "maxProbPropagationDistance must be >= 0");

        this.activityProfileStateIterator = activityProfileStateIterator;
        this.readHeader = readHeader;
        this.minRegionSize = minRegionSize;
        this.maxRegionSize = maxRegionSize;
        this.assemblyRegionPadding = assemblyRegionPadding;

        this.readyRegion = null;
        this.pendingRegions = new ArrayDeque<>();
        this.activityProfile = new BandPassActivityProfile(maxProbPropagationDistance, activeProbThreshold, BandPassActivityProfile.MAX_FILTER_SIZE, BandPassActivityProfile.DEFAULT_SIGMA, readHeader);

        readyRegion = loadNextAssemblyRegion();
    }

    @Override
    public boolean hasNext() {
        return readyRegion != null;
    }

    @Override
    public AssemblyRegion next() {
        if ( ! hasNext() ) {
            throw new NoSuchElementException("next() called when there were no more elements");
        }

        final AssemblyRegion toReturn = readyRegion;
        readyRegion = loadNextAssemblyRegion();
        return toReturn;
    }

    private AssemblyRegion loadNextAssemblyRegion() {
        AssemblyRegion nextRegion = null;

        while ( activityProfileStateIterator.hasNext() && nextRegion == null ) {
            final ActivityProfileState profile = activityProfileStateIterator.next();

            // Pop any new pending regions off of the activity profile. These pending regions will not become ready
            // until we've traversed all the reads that belong in them.
            //
            // Ordering matters here: need to check for forceConversion before adding current pileup to the activity profile
            if ( ! activityProfile.isEmpty() ) {
                final boolean forceConversion = profile.getLoc().getStart() != activityProfile.getEnd() + 1;
                pendingRegions.addAll(activityProfile.popReadyAssemblyRegions(assemblyRegionPadding, minRegionSize, maxRegionSize, forceConversion));
            }

            // Add the current pileup to the activity profile
            activityProfile.add(profile);

            // A pending region only becomes ready once our locus iterator has advanced beyond the end of its extended span
            // (this ensures that we've loaded all reads that belong in the new region)
            if ( ! pendingRegions.isEmpty() && IntervalUtils.isAfter(profile.getLoc(), pendingRegions.peek().getPaddedSpan(), readHeader.getSequenceDictionary()) ) {
                nextRegion = pendingRegions.poll();
            }
        }

        // When we run out of loci, close out the activity profile, and close out any remaining pending regions one at a time
        // It may require multiple invocations before the pendingRegions queue is cleared out.
        if ( ! activityProfileStateIterator.hasNext() ) {

            if ( ! activityProfile.isEmpty() ) {
                // Pop the activity profile a final time with forceConversion == true
                pendingRegions.addAll(activityProfile.popReadyAssemblyRegions(assemblyRegionPadding, minRegionSize, maxRegionSize, true));
            }

            // Grab the next pending region if there is one, unless we already have a region ready to go
            // (could happen if we finalize a region on the same iteration that we run out of loci in the locus iterator)
            if ( ! pendingRegions.isEmpty() && nextRegion == null ) {
                nextRegion = pendingRegions.poll();
            }
        }

        return nextRegion;
    }

    @Override
    public void remove() {
        throw new UnsupportedOperationException("remove() not supported by AssemblyRegionIterator");
    }
}
