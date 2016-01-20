package org.broadinstitute.hellbender.utils.downsampling;

import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.LinkedList;
import java.util.List;

/**
 * Pass-Through Downsampler: Implementation of the ReadsDownsampler interface that does no
 * downsampling whatsoever, and instead simply "passes-through" all the reads it's given.
 * Useful for situations where you want to disable downsampling, but still need to use
 * the downsampler interface.
 */
public final class PassThroughDownsampler<T extends GATKRead> extends ReadsDownsampler<T> {

    private LinkedList<T> selectedReads;

    public PassThroughDownsampler() {
        clearItems();
    }

    @Override
    public void submit(final T newRead ) {
        Utils.nonNull(newRead);
        // All reads pass-through, no reads get downsampled
        selectedReads.add(newRead);
    }

    @Override
    public boolean hasFinalizedItems() {
        return ! selectedReads.isEmpty();
    }

    /**
     * Note that this list is a linked list and so doesn't support fast random access
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
        selectedReads = new LinkedList<>();
    }

    @Override
    public boolean requiresCoordinateSortOrder() {
        return false;
    }

    @Override
    public void signalNoMoreReadsBefore(final T read ) {
        Utils.nonNull(read);
        // NO-OP
    }
}
