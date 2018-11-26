package org.broadinstitute.hellbender.utils.downsampling;

import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.util.*;

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
     * Random generator that should be reset for each location
     */
    private Random randomGenerator;


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
     *
     * @param useOwnRandomGenerator if true, this downsampler will rely on its own internal random
     *                              generator object that will have its seed reset. When in this mode
     *                              in order to keep the downsampling deterministic, simply supply
     *                              the read corresponding to the next start position in order to
     *                              reset the seed.
     */
    public ReservoirDownsampler(final int targetSampleSize, final boolean expectFewOverflows, final boolean useOwnRandomGenerator) {
        if ( targetSampleSize <= 0 ) {
            throw new IllegalArgumentException("Cannot do reservoir downsampling with a sample size <= 0");
        }

        this.targetSampleSize = targetSampleSize;
        this.expectFewOverflows = expectFewOverflows;
        this.randomGenerator = useOwnRandomGenerator? new Random(): null;
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
        this(targetSampleSize, false, false);
    }

    @Override
    public void submit ( final GATKRead newRead ) {
        Utils.nonNull(newRead, "newRead");

        // Only count reads that are actually eligible for discarding for the purposes of the reservoir downsampling algorithm
        totalReadsSeen++;

        if ( totalReadsSeen <= targetSampleSize ) {
            reservoir.add(newRead);
        } else {
            if ( isLinkedList ) {
                reservoir = new ArrayList<>(reservoir);
                isLinkedList = false;
            }

            final int randomSlot = (randomGenerator == null ? Utils.getRandomGenerator() : randomGenerator).nextInt(totalReadsSeen);
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

    /**
     * @param readSeed If supplied, the start position is added to the {@Link Utils#getGatkDefaultRandomSeed} to
     *                 reset the random seed for the random generator. This is used to ensure that the resivoir
     *                 downsampler is deterministic in its downsampling.
     */
    public List<GATKRead> consumeFinalizedItems(final GATKRead readSeed) {
        // Use the start position for the read as a seed for the random generator
        randomGenerator.setSeed(Utils.getGatkDefaultRandomSeed() +
                (ReadUtils.readHasNoAssignedPosition(readSeed) ? 0 : readSeed.getContig().hashCode() << 32 + readSeed.getStart()));
        return consumeFinalizedItems();
    }

    @Override
    public List<GATKRead> consumeFinalizedItems() {
        if (hasFinalizedItems()) {
            // pass reservoir by reference rather than make a copy, for speed
            final List<GATKRead> downsampledItems = reservoir;
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
        // if we aren't expecting many overflows, allocate a linked list not an arraylist
        reservoir = expectFewOverflows ? new LinkedList<>() : new ArrayList<>(targetSampleSize);

        // it's a linked list if we allocate one
        isLinkedList = expectFewOverflows;

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
