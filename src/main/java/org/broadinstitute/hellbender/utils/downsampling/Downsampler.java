/*
* Copyright 2012-2015 Broad Institute, Inc.
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.gatk.utils.downsampling;

import java.util.Collection;
import java.util.List;

/**
 * The basic downsampler API, with no reads-specific operations.
 *
 * Downsamplers that extend this class rather than the ReadsDownsampler class can handle
 * any kind of item, however they cannot be wrapped within a DownsamplingReadsIterator or a
 * PerSampleDownsamplingReadsIterator.
 *
 * @author David Roazen
 */
public abstract class Downsampler<T> {

    /**
     * Number of items discarded by this downsampler since the last call to resetStats()
     */
    protected int numDiscardedItems = 0;

    /**
     * Submit one item to the downsampler for consideration. Some downsamplers will be able to determine
     * immediately whether the item survives the downsampling process, while others will need to see
     * more items before making that determination.
     *
     * @param item the individual item to submit to the downsampler for consideration
     */
    public abstract void submit( final T item );

    /**
     * Submit a collection of items to the downsampler for consideration. Should be equivalent to calling
     * submit() on each individual item in the collection.
     *
     * @param items the collection of items to submit to the downsampler for consideration
     */
    public void submit( final Collection<T> items ) {
        if ( items == null ) {
            throw new IllegalArgumentException("submitted items must not be null");
        }

        for ( final T item : items ) {
            submit(item);
        }
    }

    /**
     * Are there items that have survived the downsampling process waiting to be retrieved?
     *
     * @return true if this downsampler has > 0 finalized items, otherwise false
     */
    public abstract boolean hasFinalizedItems();

    /**
     * Return (and *remove*) all items that have survived downsampling and are waiting to be retrieved.
     *
     * @return a list of all finalized items this downsampler contains, or an empty list if there are none
     */
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
     * Used to tell the downsampler that no more items will be submitted to it, and that it should
     * finalize any pending items.
     */
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

    /**
     * Indicates whether an item should be excluded from elimination during downsampling. By default,
     * all items representing reduced reads are excluded from downsampling, but individual downsamplers
     * may override if they are able to handle reduced reads correctly. Downsamplers should check
     * the return value of this method before discarding an item.
     *
     * @param item The item to test
     * @return true if the item should not be subject to elimination during downsampling, otherwise false
     */
    protected boolean doNotDiscardItem( final Object item ) {
        return false;
    }
}
