package org.broadinstitute.hellbender.tools.spark.transforms.markduplicates;

import com.google.api.client.util.Lists;
import com.google.common.collect.ImmutableList;
import htsjdk.samtools.*;
import org.apache.spark.SparkContext;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;
import org.broadinstitute.hellbender.engine.spark.datasources.ReadsSparkSink;
import org.broadinstitute.hellbender.utils.read.*;
import org.broadinstitute.hellbender.utils.read.markduplicates.MarkDuplicatesScoringStrategy;
import org.broadinstitute.hellbender.utils.read.markduplicates.ReadsKey;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.test.SamAssertionUtils;
import org.testng.Assert;
import org.testng.annotations.Test;
import scala.Tuple2;
import scala.collection.Seq;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

public class MarkDuplicatesSparkUtilsUnitTest extends GATKBaseTest {
    @Test(groups = "spark")
    public void testSpanningIterator() {
        check(Collections.emptyIterator(), Collections.emptyList());
        check(ImmutableList.of(pair(1, "a")).iterator(),
                ImmutableList.of(pairIterable(1, "a")));
        check(ImmutableList.of(pair(1, "a"), pair(1, "b")).iterator(),
                ImmutableList.of(pairIterable(1, "a", "b")));
        check(ImmutableList.of(pair(1, "a"), pair(2, "b")).iterator(),
                ImmutableList.of(pairIterable(1, "a"), pairIterable(2, "b")));
        check(ImmutableList.of(pair(1, "a"), pair(1, "b"), pair(2, "c")).iterator(),
                ImmutableList.of(pairIterable(1, "a", "b"), pairIterable(2, "c")));
        check(ImmutableList.of(pair(1, "a"), pair(2, "b"), pair(2, "c")).iterator(),
                ImmutableList.of(pairIterable(1, "a"), pairIterable(2, "b", "c")));
        check(ImmutableList.of(pair(1, "a"), pair(2, "b"), pair(1, "c")).iterator(),
                ImmutableList.of(pairIterable(1, "a"), pairIterable(2, "b"), pairIterable(1, "c")));
    }


    private String getReadGroupId(final SAMFileHeader header, final int index) {
        return header.getReadGroups().get(index).getReadGroupId();
    }

    private static <K, V> void check(Iterator<Tuple2<K, V>> it, List<Tuple2<K, Iterable<V>>> expected) {
        Iterator<Tuple2<K, Iterable<V>>> spanning = MarkDuplicatesSparkUtils.spanningIterator(it);
        ArrayList<Tuple2<K, Iterable<V>>> actual = Lists.newArrayList(spanning);
        Assert.assertEquals(actual, expected);
    }

    private static <K, V> Tuple2<K, V> pair(K key, V value) {
        return new Tuple2<>(key, value);
    }

    private static Tuple2<Integer, Iterable<String>> pairIterable(int i, String... s) {
        return new Tuple2<>(i, ImmutableList.copyOf(s));
    }

    private static Tuple2<String, Iterable<GATKRead>> pairIterable(String key, GATKRead... reads) {
        return new Tuple2<>(key, ImmutableList.copyOf(reads));
    }

}
