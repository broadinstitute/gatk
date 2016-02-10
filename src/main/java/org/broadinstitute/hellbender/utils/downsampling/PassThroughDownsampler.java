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

import java.util.LinkedList;
import java.util.List;

/**
 * Pass-Through Downsampler: Implementation of the ReadsDownsampler interface that does no
 * downsampling whatsoever, and instead simply "passes-through" all the reads it's given.
 * Useful for situations where you want to disable downsampling, but still need to use
 * the downsampler interface.
 *
 * @author David Roazen
 */
public class PassThroughDownsampler<T extends SAMRecord> extends ReadsDownsampler<T> {

    private LinkedList<T> selectedReads;

    public PassThroughDownsampler() {
        clearItems();
    }

    @Override
    public void submit( T newRead ) {
        // All reads pass-through, no reads get downsampled
        selectedReads.add(newRead);
    }

    @Override
    public boolean hasFinalizedItems() {
        return ! selectedReads.isEmpty();
    }

    /**
     * Note that this list is a linked list and so doesn't support fast random access
     * @return
     */
    @Override
    public List<T> consumeFinalizedItems() {
        // pass by reference rather than make a copy, for speed
        final List<T> downsampledItems = selectedReads;
        clearItems();
        return downsampledItems;
    }

    @Override
    public boolean hasPendingItems() {
        return false;
    }

    @Override
    public T peekFinalized() {
        return selectedReads.isEmpty() ? null : selectedReads.getFirst();
    }

    @Override
    public T peekPending() {
        return null;
    }

    @Override
    public int size() {
        return selectedReads.size();
    }

    @Override
    public void signalEndOfInput() {
        // NO-OP
    }

    @Override
    public void clearItems() {
        selectedReads = new LinkedList<T>();
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
