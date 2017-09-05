package org.broadinstitute.hellbender.tools.spark.sv.utils;

import htsjdk.samtools.util.Locatable;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.tools.spark.sv.utils.ShardPartitioner;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.collections.IntervalsSkipList;
import org.junit.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;
import java.util.stream.IntStream;

/**
 * Unit tests for {@link ShardPartitioner}.
 */
public final class ShardPartitionerUnitTest {

    @Test(dataProvider = "shardPartitionerData")
    public void testShardPartitioner(final int numberOfPartitions, final Collection<SimpleInterval> shards) {

        final ShardPartitioner partitioner = new ShardPartitioner(new IntervalsSkipList<>(shards), numberOfPartitions);
        Assert.assertSame(partitioner.numPartitions(), numberOfPartitions);

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
            final int[] accumulateSizes = new int[contigCount];
            for (int j = 1; j < contigCount; j++) {
                accumulateSizes[j] = accumulateSizes[j - 1] + contigSizes[j - 1];
            }
            final double depth = rdn.nextDouble() * 2;
            long totalBases = 0;
            while (totalBases < accumulateSizes[contigCount - 1] * depth) {
                int start = rdn.nextInt(accumulateSizes[contigCount - 1]);
                int contigIndex = 0;
                for (int j = 0; j < contigCount - 1; j++, contigIndex++) {
                    if (start < accumulateSizes[j + 1]) {
                        break;
                    }
                    start -= contigSizes[j];
                }
                final int length = 50;
                final int end = Math.min(start + length, contigSizes[contigIndex]);
                final SimpleInterval si = new SimpleInterval("seq" + contigIndex, start + 1, end);
                if (intervals.add(si)) {
                    totalBases += length;
                }
            }
            final int numberOfPartions = rdn.nextInt((int) Math.floor(Math.log(intervals.size()))) + 1;
            result.add(new ImmutablePair<>(numberOfPartions, intervals));
        }
        return result.stream().map(l -> new Object[] { l }).toArray(Object[][]::new);
    }


}
