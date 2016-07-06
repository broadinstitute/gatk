package org.broadinstitute.hellbender.utils.downsampling;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadCoordinateComparator;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.util.ArrayList;
import java.util.List;


/**
 * PositionalDownsampler: Downsample each stack of reads at each alignment start to a size <= a target coverage
 * using a {@link ReservoirDownsampler}. Stores only O(target coverage) reads in memory at any given time,
 * provided the client regularly calls {@link #consumeFinalizedItems}.
 *
 * Unmapped reads with assigned positions are subject to downsampling in the same way as mapped reads,
 * but unmapped reads without assigned positions are not subject to downsampling.
 */
public final class PositionalDownsampler extends ReadsDownsampler {

    private final ReservoirDownsampler reservoir;

    private final SAMFileHeader header;

    private GATKRead previousRead;

    private List<GATKRead> finalizedReads;

    /**
     * Construct a PositionalDownsampler
     *
     * @param targetCoverage Maximum number of reads that may share any given alignment start position. Must be > 0
     * @param header SAMFileHeader to use to determine contig ordering. Non-null.
     */
    public PositionalDownsampler( final int targetCoverage, final SAMFileHeader header ) {
        Utils.validateArg(targetCoverage > 0, "targetCoverage must be > 0");
        Utils.nonNull(header);

        this.reservoir = new ReservoirDownsampler(targetCoverage);
        this.finalizedReads = new ArrayList<>();
        this.header = header;
        clearItems();
        resetStats();
    }

    @Override
    public void submit( final GATKRead newRead ) {
        Utils.nonNull(newRead, "newRead");

        // If we've moved to a new position, finalize the reads currently in our reservoir.
        handlePositionalChange(newRead);

        // Pass-through reads that have no assigned position, to avoid downsampling all unmapped reads
        // to the targetCoverage. Unmapped reads that do have an assigned position *will* be subject to
        // downsampling, however.
        if ( ReadUtils.readHasNoAssignedPosition(newRead) ) {
            finalizedReads.add(newRead);
        }
        else {
            final int reservoirPreviouslyDiscardedItems = reservoir.getNumberOfDiscardedItems();
            reservoir.submit(newRead);
            incrementNumberOfDiscardedItems(reservoir.getNumberOfDiscardedItems() - reservoirPreviouslyDiscardedItems);
        }

        previousRead = newRead;
    }

    private void handlePositionalChange( final GATKRead newRead ) {
        // Use ReadCoordinateComparator to determine whether we've moved to a new start position.
        // ReadCoordinateComparator will correctly distinguish between purely unmapped reads and unmapped reads that
        // are assigned a nominal position.
        if ( previousRead != null && ReadCoordinateComparator.compareCoordinates(previousRead, newRead, header) != 0 ) {
            if ( reservoir.hasFinalizedItems() ) {
                finalizeReservoir();
            }
        }
    }

    private void finalizeReservoir() {
        finalizedReads.addAll(reservoir.consumeFinalizedItems());
        reservoir.resetStats();
    }

    @Override
    public boolean hasFinalizedItems() {
        return ! finalizedReads.isEmpty();
    }

    @Override
    public List<GATKRead> consumeFinalizedItems() {
        final List<GATKRead> toReturn = finalizedReads;
        finalizedReads = new ArrayList<>();
        return toReturn;
    }

    @Override
    public boolean hasPendingItems() {
        // The finalized items in the ReservoirDownsampler are pending items from the perspective of the
        // enclosing PositionalDownsampler
        return reservoir.hasFinalizedItems();
    }

    @Override
    public GATKRead peekFinalized() {
        return finalizedReads.isEmpty() ? null : finalizedReads.get(0);
    }

    @Override
    public GATKRead peekPending() {
        // The finalized items in the ReservoirDownsampler are pending items from the perspective of the
        // enclosing PositionalDownsampler
        return reservoir.peekFinalized();
    }

    @Override
    public int size() {
        return finalizedReads.size() + reservoir.size();
    }

    @Override
    public void signalEndOfInput() {
        finalizeReservoir();
    }

    @Override
    public void clearItems() {
        reservoir.clearItems();
        reservoir.resetStats();
        finalizedReads.clear();
        previousRead = null;
    }

    @Override
    public boolean requiresCoordinateSortOrder() {
        return true;
    }

    @Override
    public void signalNoMoreReadsBefore( final GATKRead read ) {
        handlePositionalChange(read);
    }
}
