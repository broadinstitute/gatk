package org.broadinstitute.hellbender.tools.coveragemodel;

import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.math3.util.FastMath;
import org.apache.spark.api.java.JavaPairRDD;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.factory.Nd4j;
import org.nd4j.linalg.indexing.NDArrayIndex;
import scala.Tuple2;

import javax.annotation.Nonnull;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Spark-related helper functions for {@link CoverageModelEMWorkspace}
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public final class CoverageModelSparkUtils {

    private CoverageModelSparkUtils() {}

    /**
     * Creates a uniform partitioning of an index space. If {@param minBlockSize} is larger than
     * {@param length}/{@param numBlocks}, fewer partitions will be created.
     *
     * Note: the partition number is used as the hash code.
     *
     * @param length the length of the index space
     * @param numBlocks number of partitions (= blocks)
     * @param minBlockSize minimum number of indices in each partition
     * @return a list of {@link LinearlySpacedIndexBlock}
     */
    public static List<LinearlySpacedIndexBlock> createLinearlySpacedIndexBlocks(final int length, final int numBlocks,
                                                                                 final int minBlockSize) {
        ParamUtils.isPositive(length, "The length of the linear space to be partitioned must be positive");
        ParamUtils.isPositive(numBlocks, "The number of blocks must be positive");
        ParamUtils.isPositive(minBlockSize, "Minimum block size must be positive");

        final int blockSize = FastMath.max(length / numBlocks, minBlockSize);
        final List<LinearlySpacedIndexBlock> blocks = new ArrayList<>();
        for (int begIndex = 0; begIndex < length; begIndex += blockSize) {
            blocks.add(new LinearlySpacedIndexBlock(begIndex, FastMath.min(begIndex + blockSize, length),
                    blocks.size()));
        }
        /* the last block might be smaller than minBlockSize; we merge them */
        while (blocks.size() > 1 && blocks.get(blocks.size() - 1).getNumElements() < minBlockSize) {
            final int newBegIndex = blocks.get(blocks.size() - 2).getBegIndex();
            final int newEndIndex = blocks.get(blocks.size() - 1).getEndIndex();
            blocks.remove(blocks.size() - 1);
            blocks.remove(blocks.size() - 1);
            blocks.add(new LinearlySpacedIndexBlock(newBegIndex, newEndIndex, blocks.size()));
        }
        return blocks;
    }

    /**
     * Partitions an {@link INDArray} along the 0th axis according to a given partitioning and returns
     * a list
     *
     * @param blocks partitioning
     * @param arr an array
     * @return list of partitioned arrays
     */
    public static List<Tuple2<LinearlySpacedIndexBlock, INDArray>> partitionINDArrayToList(@Nonnull final List<LinearlySpacedIndexBlock> blocks,
                                                                                           @Nonnull final INDArray arr) {
        return blocks.stream().map(block ->
                new Tuple2<>(block, arr.get(NDArrayIndex.interval(block.getBegIndex(), block.getEndIndex()))))
                .collect(Collectors.toList());
    }

    /**
     * Partitions an {@link INDArray} along the 0th axis according to a given partitioning and returns
     * a map
     *
     * @param blocks partitioning
     * @param arr an array
     * @return map of partitions to partitioned arrays
     */
    public static Map<LinearlySpacedIndexBlock, INDArray> partitionINDArrayToMap(@Nonnull final List<LinearlySpacedIndexBlock> blocks,
                                                                                 @Nonnull final INDArray arr) {
        return blocks.stream().map(block ->
                new Tuple2<>(block, arr.get(NDArrayIndex.interval(block.getBegIndex(), block.getEndIndex()))))
                .collect(Collectors.toMap(p -> p._1, p -> p._2));
    }

    /**
     * Assembles INDArray blocks in a {@code JavaPairRDD<LinearlySpacedIndexBlock, INDArray>} by collecting
     * them to a local list, sorting them based on their keys ({@link LinearlySpacedIndexBlock}), and
     * concatenating them along a given axis.
     *
     * @param blocksPairRDD an instance of {@code JavaPairRDD<LinearlySpacedIndexBlock, INDArray>}
     * @param axis axis to concat along
     * @return an instance of {@link INDArray}
     */
    public static INDArray assembleINDArrayBlocksFromRDD(@Nonnull JavaPairRDD<LinearlySpacedIndexBlock, INDArray> blocksPairRDD,
                                                         final int axis) {
        final List<INDArray> sortedBlocks = blocksPairRDD.collect().stream()
                /* sort according to begin index of each target-space partition */
                .sorted(Comparator.comparingInt(t -> t._1.getBegIndex()))
                /* remove the keys */
                .map(p -> p._2)
                /* collect to a list */
                .collect(Collectors.toList());
        /* concatenate */
        return Nd4j.concat(axis, sortedBlocks.toArray(new INDArray[sortedBlocks.size()]));
    }

    /**
     * Assemble INDArray blocks in a {@code Collection<? extends Pair<LinearlySpacedIndexBlock, INDArray>>}
     * by sorting them based on their keys ({@link LinearlySpacedIndexBlock}), and concatenating them along
     * a given axis.
     *
     * @param blocksCollection an instance of {@code Collection<? extends Pair<LinearlySpacedIndexBlock, INDArray>}
     * @param axis axis to concat along
     * @return an instance of {@link INDArray}
     */
    public static INDArray assembleINDArrayBlocksFromCollection(@Nonnull Collection<? extends Pair<LinearlySpacedIndexBlock, INDArray>> blocksCollection,
                                                                final int axis) {
        Utils.nonNull(blocksCollection, "The collection of blocks must be non-null");
        Utils.validateArg(LinearlySpacedIndexBlock.isNonOverlappingFullyCovering(blocksCollection.stream()
                .map(Pair<LinearlySpacedIndexBlock, INDArray>::getKey)
                .collect(Collectors.toList())), "Either the collection of linearly-spaced blocks is not fully covering" +
                " or contains overlapping elements");
        final List<INDArray> sortedBlocks = blocksCollection.stream()
                /* sort according to begin index of each target-space partition */
                .sorted(Comparator.comparingInt(p -> p.getKey().getBegIndex()))
                /* remove the keys */
                .map(Pair<LinearlySpacedIndexBlock, INDArray>::getValue)
                /* collect to a list */
                .collect(Collectors.toList());
        /* concatenate */
        return Nd4j.concat(axis, sortedBlocks.toArray(new INDArray[sortedBlocks.size()]));
    }

    /**
     * Takes a collection of pairs of linearly-space blocks and partially-ordered lists, sorts the collection with
     * respect to the begin index of the blocks and concatenates the partially-ordered lists to create a fully-ordered
     * list.
     *
     * Note: the blocks are assumed to be unique, non-overlapping, and fully covering
     *
     * @param blockifiedCollection a collection of pairs of {@link LinearlySpacedIndexBlock} and a list of type {@code T}
     * @param <T> type of the list
     * @return the fully ordered concatenated list
     */
    public static <T> List<T> assembleBlockifiedCollection(@Nonnull final Collection<? extends Pair<LinearlySpacedIndexBlock,
            List<T>>> blockifiedCollection) {
        Utils.nonNull(blockifiedCollection, "The collection of blocks must be non-null");
        Utils.validateArg(LinearlySpacedIndexBlock.isNonOverlappingFullyCovering(blockifiedCollection.stream()
                .map(Pair<LinearlySpacedIndexBlock, List<T>>::getKey)
                .collect(Collectors.toList())), "Either the collection of linearly-spaced blocks is not fully covering" +
                " or contains overlapping elements");
        return blockifiedCollection.stream()
                /* sort by block begin index */
                .sorted(Comparator.comparingInt(p -> p.getKey().getBegIndex()))
                /* get rid of keys */
                .map(Pair<LinearlySpacedIndexBlock, List<T>>::getValue)
                /* flat map each list and collect to a single list */
                .flatMap(List::stream).collect(Collectors.toList());
    }
}
