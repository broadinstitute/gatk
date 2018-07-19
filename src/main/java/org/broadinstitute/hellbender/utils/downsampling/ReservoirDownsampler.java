package org.broadinstitute.hellbender.utils.downsampling;

import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;

/**
 * Reservoir Downsampler: Selects n reads out of a stream whose size is not known in advance, with
 * every read in the stream having an equal chance of being selected for inclusion.
 *
 * An implementation of "Algorithm R" from the paper "Random Sampling with a Reservoir" (Jeffrey Scott Vitter, 1985)
 *
 * @author David Roazen
 */
public final class ReservoirDownsampler extends ReadsDownsampler {

    /**
     * size of our reservoir -- ie., the maximum number of reads from the stream that will be retained
     * (not including any undiscardable items)
     */
    private final int targetSampleSize;

    /**
     * At times this can be a linked list or an array list, depending on how we're accessing the
     * data and whether or not we're expecting few overflows
     */
    private List<GATKRead> reservoir;

    /**
     * Count of the number of reads seen. Used by the reservoir downsampling
     * algorithm to ensure that all reads have an equal chance of making it into the reservoir.
     */
    private int totalReadsSeen;


    /**
     * Construct a ReservoirDownsampler
     *  @param targetSampleSize Size of the reservoir used by this downsampler.
     *
     *
     */
    public ReservoirDownsampler(final int targetSampleSize) {
        if ( targetSampleSize <= 0 ) {
            throw new IllegalArgumentException("Cannot do reservoir downsampling with a sample size <= 0");
        }

        this.targetSampleSize = targetSampleSize;
        reservoir = new ArrayList<>(targetSampleSize);
        clearItems();
        resetStats();
    }

    @Override
    public void submit ( final GATKRead newRead ) {
        Utils.nonNull(newRead, "newRead");

        // Only count reads that are actually eligible for discarding for the purposes of the reservoir downsampling algorithm
        totalReadsSeen++;

        if ( totalReadsSeen <= targetSampleSize ) {
            reservoir.add(newRead);
        } else {
            final int randomSlot = Utils.getRandomGenerator().nextInt(totalReadsSeen);
            if ( randomSlot < targetSampleSize ) {
                reservoir.set(randomSlot, newRead);
            }
            incrementNumberOfDiscardedItems(1);
        }
    }

    @Override
    public boolean hasFinalizedItems() {
        return ! reservoir.isEmpty();
    }

    @Override
    public List<GATKRead> consumeFinalizedItems() {
        if (hasFinalizedItems()) {
            // one could either:
            // 1) return copy the reservoir and then clear (but not reallocate) in clearItems, or
            // 2) return the reservoir by reference and reallocate a new ArrayList in clearItems
            // 2) is worse because it always reallocates the target sample size, instead of just the reads that are present
            final List<GATKRead> downsampledItems = new ArrayList<>(reservoir);
            clearItems();
            return downsampledItems;
        } else {
            // if there's nothing here, don't bother allocating a new list
            return Collections.emptyList();
        }
    }

    @Override
    public boolean hasPendingItems() {
        return false;
    }

    @Override
    public GATKRead peekFinalized() {
        return reservoir.isEmpty() ? null : reservoir.get(0);
    }

    @Override
    public GATKRead peekPending() {
        return null;
    }

    @Override
    public int size() {
        return reservoir.size();
    }

    @Override
    public void signalEndOfInput() {
        // NO-OP
    }

    /**
     * Clear the data structures used to hold information
     */
    @Override
    public void clearItems() {
        reservoir.clear();

        // an internal stat used by the downsampling process, so not cleared by resetStats() below
        totalReadsSeen = 0;
    }

    @Override
    public boolean requiresCoordinateSortOrder() {
        return false;
    }

    @Override
    public void signalNoMoreReadsBefore(final GATKRead read ) {
        Utils.nonNull(read);
        // NO-OP
    }
}
