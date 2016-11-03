package org.broadinstitute.hellbender.utils.downsampling;

import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.iterators.PushPullTransformer;

import java.util.Collection;
import java.util.List;

/**
 * The basic downsampler API, with no reads-specific operations.
 */
public abstract class Downsampler<T> implements PushPullTransformer<T> {

    /**
     * Number of items discarded by this downsampler since the last call to resetStats()
     */
    private int numDiscardedItems = 0;

    /**
     * Submit one item to the downsampler for consideration. Some downsamplers will be able to determine
     * immediately whether the item survives the downsampling process, while others will need to see
     * more items before making that determination.
     *
     * @param item the individual item to submit to the downsampler for consideration
     */
    @Override
    public abstract void submit( final T item );

    /**
     * Submit a collection of items to the downsampler for consideration. Should be equivalent to calling
     * submit() on each individual item in the collection.
     *
     * @param items the collection of items to submit to the downsampler for consideration
     */
    @Override
    public void submit(final Collection<T> items) {
        Utils.nonNull(items, "submitted items must not be null");

        for ( final T item : items ) {
            submit(item);
        }
    }

    /**
     * Are there items that have survived the downsampling process waiting to be retrieved?
     *
     * @return true if this downsampler has > 0 finalized items, otherwise false
     */
    @Override
    public abstract boolean hasFinalizedItems();

    /**
     * Return (and *remove*) all items that have survived downsampling and are waiting to be retrieved.
     *
     * @return a list of all finalized items this downsampler contains, or an empty list if there are none
     */
    @Override
    public abstract List<T> consumeFinalizedItems();

    /**
     * Are there items stored in this downsampler that it doesn't yet know whether they will
     * ultimately survive the downsampling process?
     *
     * @return true if this downsampler has > 0 pending items, otherwise false
     */
    public abstract boolean hasPendingItems();

    /**
     * Peek at the first finalized item stored in this downsampler (or null if there are no finalized items)
     *
     * @return the first finalized item in this downsampler (the item is not removed from the downsampler by this call),
     *         or null if there are none
     */
    public abstract T peekFinalized();

    /**
     * Peek at the first pending item stored in this downsampler (or null if there are no pending items)
     *
     * @return the first pending item stored in this downsampler (the item is not removed from the downsampler by this call),
     *         or null if there are none
     */
    public abstract T peekPending();

    /**
     * Get the current number of items in this downsampler
     *
     * This should be the best estimate of the total number of elements that will come out of the downsampler
     * were consumeFinalizedItems() to be called immediately after this call.  In other words it should
     * be number of finalized items + estimate of number of pending items that will ultimately be included as well.
     *
     * @return a positive integer
     */
    public abstract int size();

    /**
     * Returns the number of items discarded (so far) during the downsampling process
     *
     * @return the number of items that have been submitted to this downsampler and discarded in the process of
     *         downsampling
     */
    public int getNumberOfDiscardedItems() {
        return numDiscardedItems;
    }

    /**
     * Increments the number of discarded items by the given value.
     *
     * Should only be called by Downsampler implementations, and so has protected access.
     *
     * @param newlyDiscardedItems amount by which to increase the number of discarded items
     */
    protected void incrementNumberOfDiscardedItems( final int newlyDiscardedItems ) {
        numDiscardedItems += newlyDiscardedItems;
    }

    /**
     * Used to tell the downsampler that no more items will be submitted to it, and that it should
     * finalize any pending items.
     */
    @Override
    public abstract void signalEndOfInput();

    /**
     * Empty the downsampler of all finalized/pending items
     */
    public abstract void clearItems();

    /**
     * Reset stats in the downsampler such as the number of discarded items *without* clearing the downsampler of items
     */
    public void resetStats() {
        numDiscardedItems = 0;
    }
}
