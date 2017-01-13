package org.broadinstitute.hellbender.tools.coveragemodel;

import org.broadinstitute.hdf5.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import javax.annotation.Nonnull;
import java.io.Serializable;
import java.util.Collection;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * This class represents a key for a contiguous block of elements in an indexed linear space.
 *
 * The begin index is include
 * The end index is exclusive
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public final class LinearSpaceBlock implements Serializable {

    private static final long serialVersionUID = 4571138983754570342L;

    /* begin index inclusive, end index exclusive */
    private final int begIndex, endIndex, numTargets;

    public LinearSpaceBlock(final int begIndex, final int endIndex) {
        this.begIndex = ParamUtils.isPositiveOrZero(begIndex, "The begin index of a block must be non-negative.");
        this.endIndex = ParamUtils.inRange(endIndex, begIndex + 1, Integer.MAX_VALUE, "The block must at least" +
                " contain one element.");
        numTargets = endIndex - begIndex;
    }

    public int getBegIndex() { return begIndex; }

    public int getEndIndex() { return endIndex; }

    public int getNumTargets() { return numTargets; }

    /**
     * Asserts that a collection of linear space blocks are non-overlapping and fully-covering (i.e. there is no
     * gap between them)
     *
     * @param blocks a collection of linear space blocks
     * @throws IllegalArgumentException if the collection is null
     * @throws AssertionError if blocks overlap or there is a gap between them
     */
    public static void assertNonOverlappingFullyCovering(@Nonnull final Collection<LinearSpaceBlock> blocks) {
        Utils.nonNull(blocks, "The collection of linear space blocks must be non-null");
        final List<LinearSpaceBlock> sortedBlocks = blocks.stream()
                .sorted(Comparator.comparingInt(LinearSpaceBlock::getBegIndex))
                .collect(Collectors.toList());
        if (!IntStream.range(0, sortedBlocks.size() - 1)
                .allMatch(idx -> sortedBlocks.get(idx).getEndIndex() == sortedBlocks.get(idx + 1).getBegIndex())) {
            throw new AssertionError("Some of the blocks in the collection either overlap or have gaps between them");
        }
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) {
            return true;
        }
        if (!(o instanceof LinearSpaceBlock)) {
            return false;
        }

        final LinearSpaceBlock block = (LinearSpaceBlock) o;
        return (begIndex == block.begIndex) && (endIndex == block.endIndex);
    }

    /**
     * The best hash code is the {@link LinearSpaceBlock#begIndex} for non-overlapping and fully-covering
     * blocks
     *
     * @return hash code
     */
    @Override
    public int hashCode() {
        return begIndex;
    }

    @Override
    public String toString() {
        return "[" + begIndex + ", " + endIndex + "]";
    }
}
