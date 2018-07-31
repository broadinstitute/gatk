package org.broadinstitute.hellbender.utils.collections;

import com.google.common.base.Stopwatch;
import com.google.common.collect.Lists;
import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.TimeUnit;

public final class IntervalsSkipListOneContigUnitTest extends GATKBaseTest {


    @DataProvider(name="intervals")
    public Object[][] intervals(){
        ArrayList<Locatable> input = Lists.newArrayList(
                new SimpleInterval("1",10,100)
        );
        ArrayList<Locatable> empty = new ArrayList<>();
        ArrayList<Locatable> manyOverlapping = Lists.newArrayList(
                new SimpleInterval("1",10,100),
                // special case: multiple intervals starting at the same place
                new SimpleInterval("1",20,50),
                new SimpleInterval("1",20,50),
                new SimpleInterval("1",20,50)
        );
        ArrayList<Locatable> mixInput = Lists.newArrayList(
                // ends before query interval
                new SimpleInterval("1",10,20),
                // ends in query interval
                new SimpleInterval("1",10,60),
                // equal to query interval
                new SimpleInterval("1",30,50),
                // covered by query interval
                new SimpleInterval("1",40,42),
                // ends after query interval
                new SimpleInterval("1",45,60),
                // starts after query interval
                new SimpleInterval("1",60,100)
        );
        ArrayList<Locatable> mixExpected = Lists.newArrayList(
                // ends in query interval
                new SimpleInterval("1",10,60),
                // equal to query interval
                new SimpleInterval("1",30,50),
                // covered by query interval
                new SimpleInterval("1",40,42),
                // ends after query interval
                new SimpleInterval("1",45,60)
        );
        // returns input single SimpleInterval, query range, expected SimpleInterval
        return new Object[][]{
                // single-point boundary cases
                new Object[]{input, new SimpleInterval("1", 10, 10), input},
                new Object[]{input, new SimpleInterval("1", 100, 100), input},
                new Object[]{input, new SimpleInterval("1", 9, 9), empty},
                new Object[]{input, new SimpleInterval("1", 11, 11), input},
                new Object[]{input, new SimpleInterval("1", 99, 99), input},
                new Object[]{input, new SimpleInterval("1", 101, 101), empty},
                // empty list boundary case
                new Object[]{empty, new SimpleInterval("1", 101, 101), empty},
                // different contig
                new Object[]{empty, new SimpleInterval("2", 101, 101), empty},
                // input exactly matches the query interval
                new Object[]{input, new SimpleInterval("1", 10, 100), input},
                // multiple intervals in the same place (potential edge case for indexing)
                new Object[]{manyOverlapping, new SimpleInterval("1", 20, 20), manyOverlapping},
                // input with multiple intervals
                new Object[]{mixInput, new SimpleInterval("1",30,50), mixExpected}

        };
    }

    @Test(dataProvider = "intervals")
    public void testOverlap(ArrayList<Locatable> input, SimpleInterval query, ArrayList<Locatable> expected) throws Exception {
        IntervalsSkipListOneContig<Locatable> ints = new IntervalsSkipListOneContig<>(input);
        List<Locatable> actual = ints.getOverlapping(query);
        Assert.assertEquals(
                actual,
                expected
        );
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testMultipleContigs() throws Exception {
        new IntervalsSkipListOneContig<>(Arrays.asList(new SimpleInterval("1",1,2), new SimpleInterval("2",1,2)));
    }

    @Test
    public void testManyIntervals() throws Exception {
        ArrayList<Locatable> si = new ArrayList<>();
        final int MAX = 10_000_000;
        for (int start = 1; start<MAX; start += 100) {
            si.add(new SimpleInterval("1",start,start+10));
            si.add(new SimpleInterval("1",start,start+200));
        }

        Stopwatch indexing = Stopwatch.createStarted();
        IntervalsSkipListOneContig<Locatable> ints = new IntervalsSkipListOneContig<>(si);
        indexing.stop();

        Stopwatch v1 = Stopwatch.createStarted();
        for (int start = 101; start<MAX; start += 5000) {
            SimpleInterval interval = new SimpleInterval("1", start + 10, start + 11);
            List<Locatable> actual = ints.getOverlappingIgnoringIndex(interval);
            Assert.assertEquals(actual.size(), 3);
            // the two that start from "start", plus the long one that starts from start-100.
            // the one that starts from start-200 ends before our test point.
            for (Locatable l : actual) {
                Assert.assertTrue(interval.overlaps(l));
            }
        }
        v1.stop();
        Stopwatch v2 = Stopwatch.createStarted();
        for (int start = 101; start<MAX; start += 5000) {
            SimpleInterval interval = new SimpleInterval("1", start + 10, start + 11);
            List<Locatable> actual = ints.getOverlapping(interval);
            Assert.assertEquals(actual.size(), 3);
            // the two that start from "start", plus the long one that starts from start-100.
            // the one that starts from start-200 ends before our test point.
            for (Locatable l : actual) {
                Assert.assertTrue(interval.overlaps(l));
            }
        }
        v2.stop();

        System.out.println("non-indexed took "+v1.elapsed(TimeUnit.MILLISECONDS)+" ms, "
                +" indexed took "+v2.elapsed(TimeUnit.MILLISECONDS)+" ms, plus "+indexing.elapsed(TimeUnit.MILLISECONDS)+" for sorting&indexing.");
    }

    @Test
    public void testLotsOfTinyIntervals() throws Exception {
        List<Locatable> input = new ArrayList<>();
        int n = 1000000;
        for (int i = 0; i < n; i++) {
            input.add(new SimpleInterval("1", 3*i+1, 3*i+2)); //1:1-2, 1:4-5, 1:7-8
        }
        final IntervalsSkipListOneContig<Locatable> skipList = new IntervalsSkipListOneContig<>(input);
        final List<Locatable> overlapping = skipList.getOverlapping(new SimpleInterval("1", 1, 3 * n + 2));
        Assert.assertEquals(input, overlapping);
    }


    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testNullCtorArg() throws Exception {
        new IntervalsSkipListOneContig<>(null);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testNullArg() throws Exception {
        List<Locatable> input = Arrays.asList(
                new SimpleInterval("1",10,100)
        );
        final IntervalsSkipListOneContig<Locatable> l = new IntervalsSkipListOneContig<>(input);
        l.getOverlapping(null);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testNotSameContig() throws Exception {
        List<Locatable> input = Arrays.asList(
                new SimpleInterval("1",10,100),
                new SimpleInterval("2",10,100)
        );
        final IntervalsSkipListOneContig<Locatable> l = new IntervalsSkipListOneContig<>(input);
    }

    @Test
    public void testQquetNotSameContig() throws Exception {
        List<Locatable> input = Arrays.asList(
                new SimpleInterval("1",10,100)
        );
        final IntervalsSkipListOneContig<Locatable> l = new IntervalsSkipListOneContig<>(input);
        final List<Locatable> res = l.getOverlappingIgnoringIndex(new SimpleInterval("2", 10, 100));
        Assert.assertEquals(res, Collections.emptyList());
    }
    @Test
    public void testEmptyInput() throws Exception {
        List<Locatable> empty = new ArrayList<>();
        final IntervalsSkipListOneContig<Locatable> l = new IntervalsSkipListOneContig<>(empty);
        Assert.assertTrue(l.getOverlapping(new SimpleInterval("", 10, 100)).isEmpty()); //try to fool it by using empty contig
        Assert.assertTrue(l.getOverlapping(new SimpleInterval("1", 10, 100)).isEmpty());
    }
}
