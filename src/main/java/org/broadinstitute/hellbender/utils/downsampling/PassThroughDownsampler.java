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
public final class PassThroughDownsampler extends ReadsDownsampler {

    private LinkedList<GATKRead> selectedReads;

    public PassThroughDownsampler() {
        clearItems();
    }

    @Override
    public void submit(final GATKRead newRead ) {
        Utils.nonNull(newRead, "newRead");
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
    public List<GATKRead> consumeFinalizedItems() {
        // pass by reference rather than make a copy, for speed
        final List<GATKRead> downsampledItems = selectedReads;
        clearItems();
        return downsampledItems;
    }

    @Override
    public boolean hasPendingItems() {
        return false;
    }

    @Override
    public GATKRead peekFinalized() {
        return selectedReads.isEmpty() ? null : selectedReads.getFirst();
    }

    @Override
    public GATKRead peekPending() {
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
    public void signalNoMoreReadsBefore(final GATKRead read ) {
        Utils.nonNull(read);
        // NO-OP
    }
}
