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
     * if true, this downsampler will be optimized for the case
     * where most of the time we won't fill up anything like the
     * targetSampleSize elements.  If this is false, we will allocate
     * internal buffers to targetSampleSize initially, which minimizes
     * the cost of allocation if we often use targetSampleSize or more
     * elements.
     */
    private final boolean expectFewOverflows;

    /**
     * At times this can be a linked list or an array list, depending on how we're accessing the
     * data and whether or not we're expecting few overflows
     */
    private List<GATKRead> reservoir;

    /**
     * Are we currently using a linked list for the reservoir?
     */
    private boolean isLinkedList;

    /**
     * Count of the number of reads seen. Used by the reservoir downsampling
     * algorithm to ensure that all reads have an equal chance of making it into the reservoir.
     */
    private int totalReadsSeen;

    /**
     * In order to guarantee that all reads have equal probability of being discarded, we need to have consumed the
     * entire input stream before any items can become finalized. All submitted items (that survive downsampling)
     * remain pending until endOfInputStream is called, at which point they become finalized.
     */
    private boolean endOfInputStream;


    /**
     * allow for deterministic behavior. used by RampedHaplotypeCaller
     */
    private boolean nonRandomReplacementMode = false;


    /**
     * Construct a ReservoirDownsampler
     *
     * @param targetSampleSize Size of the reservoir used by this downsampler.
     *
     * @param expectFewOverflows if true, this downsampler will be optimized for the case
     *                           where most of the time we won't fill up anything like the
     *                           targetSampleSize elements.  If this is false, we will allocate
     *                           internal buffers to targetSampleSize initially, which minimizes
     *                           the cost of allocation if we often use targetSampleSize or more
     *                           elements.
     */
    public ReservoirDownsampler(final int targetSampleSize, final boolean expectFewOverflows ) {
        if ( targetSampleSize <= 0 ) {
            throw new IllegalArgumentException("Cannot do reservoir downsampling with a sample size <= 0");
        }

        this.targetSampleSize = targetSampleSize;
        this.expectFewOverflows = expectFewOverflows;
        clearItems();
        resetStats();
    }

    /**
     * Construct a ReservoirDownsampler
     *
     * @param targetSampleSize Size of the reservoir used by this downsampler. Number of items retained
     *                         after downsampling will be min(totalReads, targetSampleSize)
     */
    public ReservoirDownsampler(final int targetSampleSize ) {
        this(targetSampleSize, false);
    }

    @Override
    public void submit ( final GATKRead newRead ) {
        Utils.nonNull(newRead, "newRead");
        // Once the end of the input stream has been seen, consumeFinalizedItems or clearItems must be called to
        // reset the state of the ReservoirDownsampler before more items can be submitted
        Utils.validate(! endOfInputStream, "attempt to submit read after end of input stream has been signaled");

        // Only count reads that are actually eligible for discarding for the purposes of the reservoir downsampling algorithm
        totalReadsSeen++;

        if ( totalReadsSeen <= targetSampleSize ) {
            reservoir.add(newRead);
        } else {
            if ( isLinkedList ) {
                reservoir = new ArrayList<>(reservoir);
                isLinkedList = false;
            }

            // replacing overwritten slot selection to a deterministic mechanism (optional)
            // should still provide similar distribution
            // original line is left below for reference
            final int randomSlot = !nonRandomReplacementMode
                                        ? Utils.getRandomGenerator().nextInt(totalReadsSeen)
                                        : Math.abs(newRead.getName().hashCode()) % totalReadsSeen;
            if ( randomSlot < targetSampleSize ) {
                reservoir.set(randomSlot, newRead);
            }
            incrementNumberOfDiscardedItems(1);
        }
    }

    @Override
    public boolean hasFinalizedItems() {
        // All items in the reservoir are pending until endOfInputStream is seen, at which point all items become finalized
        return endOfInputStream && !reservoir.isEmpty();
    }

    @Override
    public List<GATKRead> consumeFinalizedItems() {
        // This method clears state (including the end of input stream flag) when called after
        // end of input stream has been signaled, but has no side effects when called
        // before end of input stream has been signaled (since in that case, the downsampling
        // process is still ongoing and we shouldn't clear pending items).

        if (hasFinalizedItems()) {
            // pass reservoir by reference rather than make a copy, for speed
            final List<GATKRead> downsampledItems = reservoir;
            clearItems();
            return downsampledItems;
        } else if ( ! endOfInputStream ) {
            // Don't call clearItems() here, since endOfInputStream is false and therefore the
            // downsampling process is still ongoing. We want to preserve existing pending items,
            // and return an empty List without side effects.
            return Collections.emptyList();
        } else {
            // This is the case where endOfInputStream == true and our reservoir is empty. We return an empty
            // list, but we also call clearItems() here for consistency with the case above where we have
            // finalized items, so that in both cases we reset the endOfInputStream flag.
            clearItems();
            return Collections.emptyList();
        }
    }

    @Override
    public boolean hasPendingItems() {
        // All items in the reservoir are pending until endOfInputStream is seen, at which point all items become finalized
        return !endOfInputStream && !reservoir.isEmpty();
    }

    @Override
    public GATKRead peekFinalized() {
        return hasFinalizedItems() ? reservoir.get(0) : null;
    }

    @Override
    public GATKRead peekPending() {
        return hasPendingItems() ? reservoir.get(0) : null;
    }

    @Override
    public int size() {
        return reservoir.size();
    }

    @Override
    public void signalEndOfInput() {
        endOfInputStream = true;
    }

    /**
     * Clear the data structures used to hold information
     */
    @Override
    public void clearItems() {
        // if we aren't expecting many overflows, allocate a linked list not an arraylist
        reservoir = expectFewOverflows ? new LinkedList<>() : new ArrayList<>(targetSampleSize);

        // it's a linked list if we allocate one
        isLinkedList = expectFewOverflows;

        // an internal stat used by the downsampling process, so not cleared by resetStats() below
        totalReadsSeen = 0;

        endOfInputStream = false;
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

    public void setNonRandomReplacementMode(boolean nonRandomReplacementMode) {
        this.nonRandomReplacementMode = nonRandomReplacementMode;
    }
}
