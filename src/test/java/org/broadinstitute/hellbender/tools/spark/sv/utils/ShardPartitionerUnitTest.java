package org.broadinstitute.hellbender.tools.spark.sv.utils;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.Locatable;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.engine.Shard;
import org.broadinstitute.hellbender.engine.ShardBoundary;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.collections.IntervalsSkipList;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Unit tests for {@link ShardPartitioner}.
 */
public final class ShardPartitionerUnitTest {

    @Test(dataProvider = "shardPartitionerData")
    public void testShardPartitioner(final int numberOfPartitions, final Collection<SimpleInterval> intervals) {

        final Map<String, List<SimpleInterval>> mergedIntervals = IntervalUtils.sortAndMergeIntervals(intervals, true);
        final List<SimpleInterval> sortedIntervals = mergedIntervals
                .values().stream().flatMap(Collection::stream).collect(Collectors.toList());
        final SAMSequenceDictionary dictionary = new SAMSequenceDictionary(
                mergedIntervals.entrySet().stream()
                        .map(e -> new SAMSequenceRecord(e.getKey(), e.getValue().stream()
                                .mapToInt(Locatable::getEnd).max().getAsInt() + 10)).collect(Collectors.toList()));
        final IntervalsSkipList<ShardBoundary> shards = new IntervalsSkipList<>(sortedIntervals.stream()
                .flatMap(si -> Shard.divideIntervalIntoShards(si, 25, 0, dictionary).stream()).collect(Collectors.toList()));
        final ShardPartitioner<SimpleInterval> partitioner = new ShardPartitioner<>(SimpleInterval.class, shards, numberOfPartitions);
        Assert.assertSame(partitioner.numPartitions(), numberOfPartitions);
        for (final String contig : dictionary.getSequences().stream().map(SAMSequenceRecord::getSequenceName).collect(Collectors.toList())) {
            for (int i =1; i < dictionary.getSequence(contig).getSequenceLength(); i++) {
                final int partition = partitioner.getPartition(new SimpleInterval(contig, i, i));
                Assert.assertTrue(partition >= 0 && partition < numberOfPartitions);
                final Optional<ShardBoundary> shard = shards.getOverlapping(new SimpleInterval(contig, i, i)).stream().findFirst();
                if (shard.isPresent()) {
                    Assert.assertEquals(partition, partitioner.getPartition(new SimpleInterval(contig, i, i)));
                } else if (i > 1) {
                    Assert.assertEquals(partition, partitioner.getPartition(new SimpleInterval(contig, i - 1, i - 1)));
                } else {
                    Assert.assertEquals(partition, partitioner.getPartition(new SimpleInterval(contig, i + 1, i + 1)));
                }
            }
        }
    }

    @DataProvider(name = "shardPartitionerData")
    public Object[][] shardPartitionerData() {
        final List<Pair<Integer, Collection<SimpleInterval>>> result = new ArrayList<>();
        final Random rdn = new Random(13);
        result.add(new ImmutablePair<>(1, Collections.singleton(new SimpleInterval("seq0", 100, 200))));
        for (int i = 0; i < 1_000; i++) {
            final Set<SimpleInterval> intervals = new HashSet<>();
            final int contigCount = rdn.nextInt(3) + 1;
            final int[] contigSizes = IntStream.range(0, contigCount).map(c -> rdn.nextInt(1000) + 10).toArray();
            final int[] accumulateSizes = new int[contigCount + 1];
            for (int j = 1; j <= contigCount; j++) {
                accumulateSizes[j] = accumulateSizes[j - 1] + contigSizes[j - 1];
            }
            final double depth = rdn.nextDouble() * 2;
            long totalBases = 0;
            do {
                final int absoluteStart = rdn.nextInt(accumulateSizes[contigCount]);
                int start = absoluteStart;
                int contigIndex = 0;
                for (int j = 1; j <= contigCount; j++, contigIndex++) {
                    if (start < contigSizes[contigIndex]) {
                        break;
                    }
                    start -= contigSizes[contigIndex];
                }
                final int length = 50;
                final int end = Math.min(start + length, contigSizes[contigIndex]);
                final SimpleInterval si = new SimpleInterval("seq" + contigIndex, start + 1, end);
                if (intervals.add(si)) {
                    totalBases += length;
                }
            } while (totalBases < accumulateSizes[contigCount] * depth);
            final int numberOfPartitions = rdn.nextInt(intervals.size()) + 1;
            result.add(new ImmutablePair<>(numberOfPartitions, intervals));
        }
        return result.stream().map(p -> new Object[] { p.getLeft().intValue(), p.getRight() }).toArray(Object[][]::new);
    }


}
