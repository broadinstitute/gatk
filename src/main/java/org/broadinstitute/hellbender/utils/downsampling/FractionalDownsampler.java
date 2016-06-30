package org.broadinstitute.hellbender.utils.downsampling;

import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.ArrayList;
import java.util.List;

/**
 * Fractional Downsampler: selects a specified fraction of the reads for inclusion.
 *
 * Since the selection is done randomly, the actual fraction of reads retained may be slightly
 * more or less than the requested fraction, depending on the total number of reads submitted.
 *
 * @author David Roazen
 */
public final class FractionalDownsampler extends ReadsDownsampler {

    private List<GATKRead> selectedReads;

    private final int cutoffForInclusion;

    private static final int RANDOM_POOL_SIZE = 10000;

    /**
     * Construct a FractionalDownsampler
     *
     * @param fraction Fraction of reads to preserve, between 0.0 (inclusive) and 1.0 (inclusive).
     *                 Actual number of reads preserved may differ randomly.
     */
    public FractionalDownsampler( final double fraction ) {
        if ( fraction < 0.0 || fraction > 1.0 ) {
            throw new IllegalArgumentException("Fraction of reads to include must be between 0.0 and 1.0, inclusive");
        }

        cutoffForInclusion = (int)(fraction * RANDOM_POOL_SIZE);
        clearItems();
        resetStats();
    }

    @Override
    public void submit( final GATKRead newRead ) {
        Utils.nonNull(newRead, "newRead");

        if ( Utils.getRandomGenerator().nextInt(10000) < cutoffForInclusion) {
            selectedReads.add(newRead);
        } else {
            incrementNumberOfDiscardedItems(1);
        }
    }

    @Override
    public boolean hasFinalizedItems() {
        return !selectedReads.isEmpty();
    }

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
        return selectedReads.isEmpty() ? null : selectedReads.get(0);
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
        selectedReads = new ArrayList<>();
    }

    @Override
    public boolean requiresCoordinateSortOrder() {
        return false;
    }

    @Override
    public void signalNoMoreReadsBefore( final GATKRead read ) {
        Utils.nonNull(read);
        // NO-OP
    }
}
