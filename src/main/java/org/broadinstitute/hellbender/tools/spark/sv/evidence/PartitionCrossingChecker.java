package org.broadinstitute.hellbender.tools.spark.sv.evidence;

import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;

/** It allows you to ask whether a given interval is near the beginning or end of the partition.
 *  This is useful in handling breakpoints for which the evidence is split across two partitions. */
public class PartitionCrossingChecker {
    private final int beginningContigID;
    private final int beginningPosition;
    private final int endingContigID;
    private final int endingPosition;

    /** A constructor that will always report onBoundary to be false. */
    public PartitionCrossingChecker() {
        beginningContigID = -1;
        beginningPosition = 0;
        endingContigID = -1;
        endingPosition = 0;
    }

    /** A constructor from the metadata with a specified boundary width. */
    public PartitionCrossingChecker( final int partitionIdx,
                                     final ReadMetadata readMetadata,
                                     final int boundaryWidth ) {
        if ( partitionIdx == 0 ) {
            beginningContigID = ReadMetadata.PartitionBounds.UNMAPPED;
            beginningPosition = -1;
        } else {
            final ReadMetadata.PartitionBounds bounds = readMetadata.getPartitionBounds(partitionIdx - 1);
            beginningContigID = bounds.getLastContigID();
            if ( beginningContigID == ReadMetadata.PartitionBounds.UNMAPPED ) {
                beginningPosition = -1;
            } else {
                beginningPosition = bounds.getLastStart() + boundaryWidth;
            }
        }

        if ( partitionIdx == readMetadata.getNPartitions() - 1 ) {
            endingContigID = ReadMetadata.PartitionBounds.UNMAPPED;
            endingPosition = Integer.MAX_VALUE;
        } else {
            final ReadMetadata.PartitionBounds bounds = readMetadata.getPartitionBounds(partitionIdx + 1);
            endingContigID = bounds.getFirstContigID();
            if ( endingContigID == ReadMetadata.PartitionBounds.UNMAPPED ) {
                endingPosition = Integer.MAX_VALUE;
            } else {
                endingPosition = bounds.getFirstStart() - boundaryWidth;
            }
        }
    }

    public boolean onBoundary( final SVInterval interval ) {
        return onLeadingBoundary(interval) || onTrailingBoundary(interval);
    }

    public boolean onLeadingBoundary( final SVInterval interval ) {
        return interval.getContig() == beginningContigID && interval.getStart() < beginningPosition;
    }

    public boolean onTrailingBoundary( final SVInterval interval ) {
        return interval.getContig() == endingContigID && interval.getEnd() >= endingPosition;
    }
}
