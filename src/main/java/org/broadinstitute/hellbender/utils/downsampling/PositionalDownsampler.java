package org.broadinstitute.hellbender.utils.downsampling;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadCoordinateComparator;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;


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

    private List<GATKRead> finalizedReads;

    private GATKRead firstReadInStride;

    private final int alignmentStartStride;

    /**
     * Construct a PositionalDownsampler
     *
     * @param targetCoverage Maximum number of reads that may share any given alignment start position. Must be > 0
     * @param header SAMFileHeader to use to determine contig ordering. Non-null.
     * @param downsampleByMappingQuality    If true, bias downsampling toward reads with higher mapping quality.
     */
    public PositionalDownsampler( final int targetCoverage, final int alignmentStartStride, final SAMFileHeader header, final boolean downsampleByMappingQuality, final int depthToIgnoreLocus) {
        Utils.validateArg(targetCoverage > 0, "targetCoverage must be > 0");
        Utils.nonNull(header);

        this.reservoir = new ReservoirDownsampler(targetCoverage, false, downsampleByMappingQuality, depthToIgnoreLocus);
        this.finalizedReads = new ArrayList<>();
        this.header = header;
        this.alignmentStartStride = alignmentStartStride;
        clearItems();
        resetStats();
    }

    /**
     * Construct a PositionalDownsampler
     *
     * @param targetCoverage Maximum number of reads that may share any given alignment start position. Must be > 0
     * @param header SAMFileHeader to use to determine contig ordering. Non-null.
     */
    public PositionalDownsampler( final int targetCoverage, final SAMFileHeader header ) {
        this(targetCoverage, 1, header, false, Integer.MAX_VALUE);
    }

    @Override
    public void submit( final GATKRead newRead ) {
        Utils.nonNull(newRead, "newRead");

        // Pass-through reads that have no assigned position, to avoid downsampling all unmapped reads
        // to the targetCoverage. Unmapped reads that do have an assigned position *will* be subject to
        // downsampling, however.
        if ( ReadUtils.readHasNoAssignedPosition(newRead) ) {
            finalizedReads.add(newRead);
        } else {
            // If we've moved to a new position, finalize the reads currently in our reservoir.
            handlePositionalChange(newRead);
            final int reservoirPreviouslyDiscardedItems = reservoir.getNumberOfDiscardedItems();
            reservoir.submit(newRead);
            incrementNumberOfDiscardedItems(reservoir.getNumberOfDiscardedItems() - reservoirPreviouslyDiscardedItems);
        }
    }

    private void handlePositionalChange( final GATKRead newRead ) {
        if (ReadUtils.readHasNoAssignedPosition(newRead)) {
            return;
        } else if (firstReadInStride == null) {
            firstReadInStride = newRead;
        } else if (!newRead.getAssignedContig().equals(firstReadInStride.getAssignedContig())
                || newRead.getAssignedStart() >= firstReadInStride.getAssignedStart() + alignmentStartStride) {
            finalizeReservoir();
            firstReadInStride = newRead;
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

        // if the stride length is > 1, different alignment start positions may be mixed in the reservoir
        return alignmentStartStride == 1 ? toReturn : toReturn.stream()
                .sorted(Comparator.comparingInt(GATKRead::getAssignedStart))
                .collect(Collectors.toList());
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
        firstReadInStride = null;
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
