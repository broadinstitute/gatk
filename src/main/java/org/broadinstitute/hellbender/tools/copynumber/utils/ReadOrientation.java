package org.broadinstitute.hellbender.tools.copynumber.utils;

import com.google.common.annotations.VisibleForTesting;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.function.Function;

/**
 * Helper class to calculate fragment center of a properly paired read
 */
@VisibleForTesting
public enum ReadOrientation {
    /**
     * Read was located on forward strand
     */
    FORWARD(read -> read.getUnclippedStart() + read.getFragmentLength() / 2,
            read -> new ImmutablePair<>(read.getUnclippedStart(), read.getUnclippedStart() + read.getFragmentLength())), //TODO check if you need to add/subtract 1
    /**
     * Read was located on reverse strand
     */
    REVERSE(read -> read.getUnclippedStart() + (read.getLength() - 1)  + read.getFragmentLength() / 2,
            read -> new ImmutablePair<>(read.getUnclippedEnd() + read.getFragmentLength(), read.getUnclippedEnd())); //TODO check if you need to add/subtract 1

    private final Function<GATKRead, Integer> readToFragmentCenterMapper;
    private final Function<GATKRead, Pair<Integer, Integer>> readToStartAndEndFragmentPositionsMapper;

    ReadOrientation(final Function<GATKRead, Integer> readToCenterMapper,
                    final Function<GATKRead, Pair<Integer, Integer>> readToStartAndEndFragmentPositionsMapper) {
        this.readToFragmentCenterMapper = readToCenterMapper;
        this.readToStartAndEndFragmentPositionsMapper = readToStartAndEndFragmentPositionsMapper;
    }

    /**
     * Get a function that maps the read to the center of the fragment.
     */
    public Function<GATKRead, Integer> getReadToFragmentCenterMapper() {
        return readToFragmentCenterMapper;
    }

    /**
     * Get a function that maps the read to a pair of integers corresponding to read's fragment start and end positions
     */
    public Function<GATKRead, Pair<Integer, Integer>> getReadToStartAndEndFragmentPositionsMapper() {
        return readToStartAndEndFragmentPositionsMapper;
    }

    /**
     * Get {@link ReadOrientation} instance corresponding to the orientation of the read.
     */
    public static ReadOrientation getReadOrientation(final GATKRead read) {
        return read.getFragmentLength() > 0 ? FORWARD : REVERSE;
    }

    /**
     * Compute center of the fragment corresponding to the read.
     */
    public static SimpleInterval getFragmentCenter(final GATKRead read) {
        final int fragmentCenter = getReadOrientation(read).getReadToFragmentCenterMapper().apply(read);
        return new SimpleInterval(read.getContig(), fragmentCenter, fragmentCenter);
    }
}
