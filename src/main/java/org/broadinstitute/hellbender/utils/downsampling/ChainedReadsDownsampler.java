package org.broadinstitute.hellbender.utils.downsampling;

import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.List;
import java.util.Optional;

/**
 * Runs two reads downsamplers in sequence: every read finalized by the first is submitted to the second, and the
 * second's finalized reads are this downsampler's output. Discard counts are summed.
 */
public final class ChainedReadsDownsampler extends ReadsDownsampler {

    private final ReadsDownsampler first;
    private final ReadsDownsampler second;

    public ChainedReadsDownsampler(final ReadsDownsampler first, final ReadsDownsampler second) {
        this.first = Utils.nonNull(first);
        this.second = Utils.nonNull(second);
    }

    private void drainFirstIntoSecond() {
        if (first.hasFinalizedItems()) {
            second.submit(first.consumeFinalizedItems());
        }
    }

    @Override
    public void submit(final GATKRead item) {
        first.submit(item);
        drainFirstIntoSecond();
    }

    @Override
    public boolean hasFinalizedItems() {
        return second.hasFinalizedItems();
    }

    @Override
    public List<GATKRead> consumeFinalizedItems() {
        return second.consumeFinalizedItems();
    }

    @Override
    public boolean hasPendingItems() {
        return first.hasPendingItems() || first.hasFinalizedItems() || second.hasPendingItems();
    }

    @Override
    public GATKRead peekFinalized() {
        return second.peekFinalized();
    }

    @Override
    public GATKRead peekPending() {
        // The pending item nearest to being emitted is the one furthest down the chain: waiting in the second
        // stage, then finalized by the first stage but not yet handed over, then still pending in the first stage.
        return Optional.ofNullable(second.peekPending())
                .or(() -> Optional.ofNullable(first.peekFinalized()))
                .orElse(first.peekPending());
    }

    @Override
    public int size() {
        return first.size() + second.size();
    }

    @Override
    public int getNumberOfDiscardedItems() {
        return first.getNumberOfDiscardedItems() + second.getNumberOfDiscardedItems();
    }

    @Override
    public void signalEndOfInput() {
        first.signalEndOfInput();
        drainFirstIntoSecond();
        second.signalEndOfInput();
    }

    @Override
    public void clearItems() {
        first.clearItems();
        second.clearItems();
    }

    @Override
    public void resetStats() {
        first.resetStats();
        second.resetStats();
    }

    @Override
    public boolean requiresCoordinateSortOrder() {
        return first.requiresCoordinateSortOrder() || second.requiresCoordinateSortOrder();
    }

    @Override
    public void signalNoMoreReadsBefore(final GATKRead read) {
        first.signalNoMoreReadsBefore(read);
        drainFirstIntoSecond();
        second.signalNoMoreReadsBefore(read);
    }
}
