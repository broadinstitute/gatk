package org.broadinstitute.hellbender.utils.iterators;


import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;

public class IntervalLocusIteratorUnitTest extends GATKBaseTest {

    @Test
    public void testSimple() {
        final SimpleInterval record_1_1_100 = new SimpleInterval("1:1-100");
        final SimpleInterval record_1_500_600 = new SimpleInterval("1:500-600");
        final SimpleInterval record_1_800_800 = new SimpleInterval("1:800-800");
        List<SimpleInterval> intervals = new ArrayList<>();
        intervals.add(record_1_1_100);
        intervals.add(record_1_500_600);
        intervals.add(record_1_800_800);

        final IntervalLocusIterator intervalLocusIterator = new IntervalLocusIterator(intervals.iterator());

        final List<SimpleInterval> newIntervals = StreamSupport.stream(Spliterators.spliteratorUnknownSize(intervalLocusIterator, Spliterator.ORDERED),
                false).collect(Collectors.toList());

        Assert.assertEquals(newIntervals.size(), 202);
        Assert.assertEquals(newIntervals.get(0).getStart(), 1);
        Assert.assertEquals(newIntervals.get(1).getStart(), 2);
        Assert.assertEquals(newIntervals.get(99).getStart(), 100);
        Assert.assertEquals(newIntervals.get(100).getStart(), 500);
        Assert.assertEquals(newIntervals.get(200).getStart(), 600);
        Assert.assertEquals(newIntervals.get(201).getStart(), 800);

        // Make sure that start = end for all intervals coming out of the LocusIntervalIterator, since each is only 1 bp
        Assert.assertTrue(newIntervals.stream().allMatch(i -> i.getEnd() == i.getStart()));
    }

    @Test
    public void testSimpleSinglePoints() {
        final SimpleInterval record_1_1_1 = new SimpleInterval("1:1-1");
        final SimpleInterval record_1_5_5 = new SimpleInterval("1:5-5");
        List<SimpleInterval> intervals = new ArrayList<>();
        intervals.add(record_1_1_1);
        intervals.add(record_1_5_5);

        final IntervalLocusIterator intervalLocusIterator = new IntervalLocusIterator(intervals.iterator());

        final List<SimpleInterval> newIntervals = StreamSupport.stream(Spliterators.spliteratorUnknownSize(intervalLocusIterator, Spliterator.ORDERED),
                false).collect(Collectors.toList());

        Assert.assertEquals(newIntervals.size(), 2);
        Assert.assertEquals(newIntervals.get(0).getStart(), 1);
        Assert.assertEquals(newIntervals.get(1).getStart(), 5);
        Assert.assertEquals(newIntervals.get(0).getEnd(), 1);
        Assert.assertEquals(newIntervals.get(1).getEnd(), 5);
    }

    @Test
    public void testEmpty() {
        final List<SimpleInterval> intervals = new ArrayList<>();

        final IntervalLocusIterator intervalLocusIterator = new IntervalLocusIterator(intervals.iterator());

        final List<SimpleInterval> newIntervals = StreamSupport.stream(Spliterators.spliteratorUnknownSize(intervalLocusIterator, Spliterator.ORDERED),
                false).collect(Collectors.toList());

        Assert.assertEquals(newIntervals.size(), 0);
    }

    @Test
    public void testSlightChallenge() {

        final SimpleInterval record_20_9999600_9999800 = new SimpleInterval("20:9999600-9999600");
        final SimpleInterval record_20_9999800_10000000 = new SimpleInterval("20:9999800-10000000");
        final SimpleInterval record_20_10000050_10000100 = new SimpleInterval("20:10000050-10000100");
        List<SimpleInterval> intervals = new ArrayList<>();
        intervals.add(record_20_9999600_9999800);
        intervals.add(record_20_9999800_10000000);
        intervals.add(record_20_10000050_10000100);

        final IntervalLocusIterator intervalLocusIterator = new IntervalLocusIterator(intervals.iterator());

        final List<SimpleInterval> newIntervals = StreamSupport.stream(Spliterators.spliteratorUnknownSize(intervalLocusIterator, Spliterator.ORDERED),
                false).collect(Collectors.toList());

        Assert.assertEquals(newIntervals.size(), 253);
    }

    @Test
    public void testSingleInterval() {
        final SimpleInterval record_20_10000098_10000101 = new SimpleInterval("20:10000098-10000101");
        final List<SimpleInterval> intervals = new ArrayList<>();
        intervals.add(record_20_10000098_10000101);
        final IntervalLocusIterator intervalLocusIterator = new IntervalLocusIterator(intervals.iterator());

        final List<SimpleInterval> newIntervals = StreamSupport.stream(Spliterators.spliteratorUnknownSize(intervalLocusIterator, Spliterator.ORDERED),
                false).collect(Collectors.toList());

        Assert.assertEquals(newIntervals.size(), 4);
    }

}
