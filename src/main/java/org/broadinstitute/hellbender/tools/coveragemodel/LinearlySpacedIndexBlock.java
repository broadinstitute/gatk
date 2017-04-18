package org.broadinstitute.hellbender.tools.coveragemodel;

import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import javax.annotation.Nonnull;
import java.io.Serializable;
import java.util.Collection;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * This class represents a contiguous block of elements in a linearly-indexed space. It is used, for instance, as a
 * lightweight "key" for pulling a contiguous block of a partitioned vector, matrix, or tensor.
 *
 * Notes:
 *
 * - The begin index is inclusive, and the end index is exclusive
 * - The hash code is based on the {@link #begIndex}
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
final class LinearlySpacedIndexBlock implements Serializable {

    private static final long serialVersionUID = 4571138983754570342L;

    private final int begIndex;
    private final int endIndex;
    private final int numElements;
    private final int hashCode;

    /**
     * Public constructor.
     *
     * @param begIndex begin index (non-negative integer)
     * @param endIndex end index (non-negative integer > begIndex)
     * @param hashCode a custom hash code (it is the user's responsibility to choose a good hash function)
     */
    public LinearlySpacedIndexBlock(final int begIndex, final int endIndex, final int hashCode) {
        Utils.validateArg(endIndex > begIndex, "The block must contain at least one element");
        this.begIndex = ParamUtils.isPositiveOrZero(begIndex, "The begin index of a block must be non-negative");
        this.endIndex = ParamUtils.isPositiveOrZero(endIndex, "The end index of a block must be non-negative");
        this.hashCode = hashCode;
        numElements = endIndex - begIndex;
    }

    public int getBegIndex() { return begIndex; }

    public int getEndIndex() { return endIndex; }

    public int getNumElements() { return numElements; }

    /**
     * Checks whether a collection of linearly spaced index blocks are non-overlapping and fully-covering (i.e. there
     * is no gap between them)
     *
     * @param blocks a collection of linear space blocks
     * @throws IllegalArgumentException if the collection is null
     */
    public static boolean isNonOverlappingFullyCovering(@Nonnull final Collection<LinearlySpacedIndexBlock> blocks) {
        Utils.nonNull(blocks, "The collection of linear space blocks must be non-null");
        final List<LinearlySpacedIndexBlock> sortedBlocks = blocks.stream()
                .sorted(Comparator.comparingInt(LinearlySpacedIndexBlock::getBegIndex))
                .collect(Collectors.toList());
        return IntStream.range(0, sortedBlocks.size() - 1)
                .allMatch(idx -> sortedBlocks.get(idx).getEndIndex() == sortedBlocks.get(idx + 1).getBegIndex());
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) {
            return true;
        }
        if (!(o instanceof LinearlySpacedIndexBlock)) {
            return false;
        }

        final LinearlySpacedIndexBlock block = (LinearlySpacedIndexBlock) o;
        return (begIndex == block.begIndex) && (endIndex == block.endIndex);
    }

    @Override
    public int hashCode() {
        return hashCode;
    }

    @Override
    public String toString() {
        return "LinearlySpacedIndexBlock{" +
                "begIndex=" + begIndex +
                ", endIndex=" + endIndex +
                ", numElements=" + numElements +
                ", hashCode=" + hashCode +
                '}';
    }
}
