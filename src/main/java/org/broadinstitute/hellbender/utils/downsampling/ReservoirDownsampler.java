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

import htsjdk.samtools.SAMRecord;
import org.broadinstitute.gatk.utils.Utils;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;

import java.util.*;

/**
 * Reservoir Downsampler: Selects n reads out of a stream whose size is not known in advance, with
 * every read in the stream having an equal chance of being selected for inclusion.
 *
 * An implementation of "Algorithm R" from the paper "Random Sampling with a Reservoir" (Jeffrey Scott Vitter, 1985)
 *
 * @author David Roazen
 */
public class ReservoirDownsampler<T extends SAMRecord> extends ReadsDownsampler<T> {

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
    private List<T> reservoir;

    /**
     * Certain items (eg., reduced reads) cannot be discarded at all during downsampling. We store
     * these items separately so as not to impact the fair selection of items for inclusion in the
     * reservoir. These items are returned (and cleared) along with any items in the reservoir in
     * calls to consumeFinalizedItems().
     */
    private List<T> undiscardableItems;

    /**
     * Are we currently using a linked list for the reservoir?
     */
    private boolean isLinkedList;

    /**
     * Count of the number of reads seen that were actually eligible for discarding. Used by the reservoir downsampling
     * algorithm to ensure that all discardable reads have an equal chance of making it into the reservoir.
     */
    private int totalDiscardableReadsSeen;


    /**
     * Construct a ReservoirDownsampler
     *
     * @param targetSampleSize Size of the reservoir used by this downsampler. Number of items retained
     *                         after downsampling will be min(totalDiscardableReads, targetSampleSize) + any
     *                         undiscardable reads (eg., reduced reads).
     *
     * @param expectFewOverflows if true, this downsampler will be optimized for the case
     *                           where most of the time we won't fill up anything like the
     *                           targetSampleSize elements.  If this is false, we will allocate
     *                           internal buffers to targetSampleSize initially, which minimizes
     *                           the cost of allocation if we often use targetSampleSize or more
     *                           elements.
     */
    public ReservoirDownsampler ( final int targetSampleSize, final boolean expectFewOverflows ) {
        if ( targetSampleSize <= 0 ) {
            throw new ReviewedGATKException("Cannot do reservoir downsampling with a sample size <= 0");
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
    public ReservoirDownsampler ( final int targetSampleSize ) {
        this(targetSampleSize, false);
    }

    @Override
    public void submit ( final T newRead ) {
        if ( doNotDiscardItem(newRead) ) {
            undiscardableItems.add(newRead);
            return;
        }

        // Only count reads that are actually eligible for discarding for the purposes of the reservoir downsampling algorithm
        totalDiscardableReadsSeen++;

        if ( totalDiscardableReadsSeen <= targetSampleSize ) {
            reservoir.add(newRead);
        }
        else {
            if ( isLinkedList ) {
                reservoir = new ArrayList<T>(reservoir);
                isLinkedList = false;
            }

            final int randomSlot = Utils.getRandomGenerator().nextInt(totalDiscardableReadsSeen);
            if ( randomSlot < targetSampleSize ) {
                reservoir.set(randomSlot, newRead);
            }
            numDiscardedItems++;
        }
    }

    @Override
    public boolean hasFinalizedItems() {
        return ! reservoir.isEmpty() || ! undiscardableItems.isEmpty();
    }

    @Override
    public List<T> consumeFinalizedItems() {
        if ( ! hasFinalizedItems() ) {
            // if there's nothing here, don't bother allocating a new list
            return Collections.emptyList();
        } else {
            // pass reservoir by reference rather than make a copy, for speed
            final List<T> downsampledItems = reservoir;
            downsampledItems.addAll(undiscardableItems);
            clearItems();
            return downsampledItems;
        }
    }

    @Override
    public boolean hasPendingItems() {
        return false;
    }

    @Override
    public T peekFinalized() {
        return ! reservoir.isEmpty() ? reservoir.get(0) : (! undiscardableItems.isEmpty() ? undiscardableItems.get(0) : null);
    }

    @Override
    public T peekPending() {
        return null;
    }

    @Override
    public int size() {
        return reservoir.size() + undiscardableItems.size();
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
        reservoir = expectFewOverflows ? new LinkedList<T>() : new ArrayList<T>(targetSampleSize);

        // there's no possibility of overflow with the undiscardable items, so we always use a linked list for them
        undiscardableItems = new LinkedList<>();

        // it's a linked list if we allocate one
        isLinkedList = expectFewOverflows;

        // an internal stat used by the downsampling process, so not cleared by resetStats() below
        totalDiscardableReadsSeen = 0;
    }

    @Override
    public boolean requiresCoordinateSortOrder() {
        return false;
    }

    @Override
    public void signalNoMoreReadsBefore( T read ) {
        // NO-OP
    }
}
