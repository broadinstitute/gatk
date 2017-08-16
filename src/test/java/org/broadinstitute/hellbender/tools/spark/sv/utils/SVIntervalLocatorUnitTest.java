package org.broadinstitute.hellbender.tools.spark.sv.utils;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.apache.commons.lang3.tuple.Triple;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;

/**
 * Unit tests for {@link SVIntervalLocator}
 */
public class SVIntervalLocatorUnitTest {

    private static final SAMSequenceDictionary TEST_DICTIONARY =
            new SAMSequenceDictionary(Arrays.asList(
                    new SAMSequenceRecord("seq1", 100),
                    new SAMSequenceRecord("seq2", 1000),
                    new SAMSequenceRecord("seq3", 10000)
            ));

    private static final SAMSequenceDictionary TEST_DICTIONARY_FOR_TREE =
            new SAMSequenceDictionary(Arrays.asList(
                    new SAMSequenceRecord("seq1", 10000),
                    new SAMSequenceRecord("seq2", 10000),
                    new SAMSequenceRecord("seq3", 10000)
            ));

    private static final SVIntervalLocator TEST_LOCATOR = SVIntervalLocator.of(TEST_DICTIONARY);

    @Test(dataProvider = "simpleToSVIntervalData")
    public void testToSVInterval(final SimpleInterval si, final SVInterval sv, final int padding) {
        Assert.assertEquals(TEST_LOCATOR.toSVInterval(si, padding), sv);
        if (padding == 0) {
            Assert.assertEquals(TEST_LOCATOR.toSVInterval(si), sv);
        }
    }

    @Test(dataProvider = "svToSimpleIntervalData")
    public void testToLocatable(final SVInterval sv, final SimpleInterval si, final int padding) {
        Assert.assertEquals(TEST_LOCATOR.toSimpleInterval(sv, padding), si);
        Assert.assertEquals(TEST_LOCATOR.toLocatable(sv, padding, SimpleInterval::new), si);
        if (padding == 0) {
            Assert.assertEquals(TEST_LOCATOR.toSimpleInterval(sv), si);
            Assert.assertEquals(TEST_LOCATOR.toLocatable(sv, SimpleInterval::new), si);
        }
    }

    @Test
    public void testAsDictionary() {
        final SAMSequenceDictionary dictionary = TEST_LOCATOR.asDictionary();
        Assert.assertEquals(dictionary.size(), TEST_DICTIONARY.size());
        for (int i = 0; i < dictionary.size(); i++) {
            Assert.assertEquals(dictionary.getSequence(i).getSequenceName(), TEST_DICTIONARY.getSequence(i).getSequenceName());
            Assert.assertEquals(dictionary.getSequence(i).getSequenceLength(), TEST_DICTIONARY.getSequence(i).getSequenceLength());
        }
    }

    @Test
    public void testToTree() {
        final SVIntervalLocator locator = SVIntervalLocator.of(TEST_DICTIONARY_FOR_TREE);
        final Random rdn = new Random(3131);
        final Set<SimpleInterval> intervals = new LinkedHashSet<>(100);
        for (int i = 0; i < 100;) {
            if (intervals.add(randomInterval(rdn, TEST_DICTIONARY_FOR_TREE, null, 1, 10000, 1, 250))) {
                i++;
            }
        }
        final SVIntervalTree<SimpleInterval> tree = intervals.stream()
                .collect(locator.toSVIntervalTree(si -> si));
        Assert.assertNotNull(tree);
        Assert.assertEquals(tree.size(), intervals.size());
        for (final SimpleInterval interval : intervals) {
            final SVInterval svInterval = locator.toSVInterval(interval);
            final SVIntervalTree.Entry<SimpleInterval> entry = tree.find(svInterval);
            Assert.assertNotNull(entry);
            Assert.assertEquals(entry.getInterval(), svInterval);
            Assert.assertEquals(entry.getValue(), interval);
        }
    }

    @DataProvider(name="simpleToSVIntervalData")
    public Object[][] simpleToSVIntervalData() {
        final Random rdn = new Random(123);
        final List<Triple<SimpleInterval, SVInterval, Integer>> result = new ArrayList<>();

        for (int i = 0; i < 100; i++) {
            final SimpleInterval interval = randomInterval(rdn, TEST_DICTIONARY, "seq2", 1, 990, 1, 100);
            final SVInterval svInterval = new SVInterval(1, interval.getStart(), interval.getEnd() + 1);
            result.add(Triple.of(interval, svInterval, 0));
            final SVInterval svInterval10Padding = new SVInterval(1,
                    Math.max(1, svInterval.getStart() - 55), Math.min(1001, svInterval.getEnd() + 55));
            result.add(Triple.of(interval, svInterval10Padding, 55));
        }
        for (int i = 0; i < 100; i++) {
            final SimpleInterval interval = randomInterval(rdn, TEST_DICTIONARY, "seq1", 1, 90, 1, 10);
            final SVInterval svInterval = new SVInterval(0, interval.getStart(), interval.getEnd() + 1);
            result.add(Triple.of(interval, svInterval, 0));
            final SVInterval svInterval10Padding = new SVInterval(0,
                    Math.max(1, svInterval.getStart() - 10), Math.min(101, svInterval.getEnd() + 10));
            result.add(Triple.of(interval, svInterval10Padding, 10));
        }
        result.add(Triple.of(new SimpleInterval("seq1", 1,1), new SVInterval(0, 1,2), 0));
        result.add(Triple.of(new SimpleInterval("seq1", 100, 100), new SVInterval(0, 100, 101), 0));
        result.add(Triple.of(new SimpleInterval("seq1", 100, 100), new SVInterval(0, 1, 101), 2000));

        result.add(Triple.of(new SimpleInterval("seq1", 10, 30), new SVInterval(0, 1, 51), 20));
        result.add(Triple.of(new SimpleInterval("seq3", 1, 10000), new SVInterval(2, 1, 10001), 0));
        result.add(Triple.of(new SimpleInterval("seq3", 1, 10000), new SVInterval(2, 1, 10001), 1));
        result.add(Triple.of(new SimpleInterval("seq3", 1, 10000), new SVInterval(2, 1, 10001), 1111));

        return result.stream()
                .map(t -> new Object[] {t.getLeft(), t.getMiddle(), t.getRight()})
                .toArray(Object[][]::new);
    }


    @DataProvider(name="svToSimpleIntervalData")
    public Object[][] svToSimpleIntervalData() {
        final Random rdn = new Random(321);
        final List<Triple<SVInterval, SimpleInterval, Integer>> result = new ArrayList<>();

        for (int i = 0; i < 100; i++) {
            final SimpleInterval interval = randomInterval(rdn, TEST_DICTIONARY, "seq2", 1, 990, 1, 100);
            final SVInterval svInterval = new SVInterval(1, interval.getStart(), interval.getEnd() + 1);
            result.add(Triple.of(svInterval, interval, 0));
            final SimpleInterval simpleInterval55Padding = new SimpleInterval("seq2",
                    Math.max(1, interval.getStart() - 55), Math.min(1000, interval.getEnd() + 55));
            result.add(Triple.of(svInterval, simpleInterval55Padding, 55));
        }
        for (int i = 0; i < 100; i++) {
            final SimpleInterval interval = randomInterval(rdn, TEST_DICTIONARY, "seq1", 1, 90, 1, 10);
            final SVInterval svInterval = new SVInterval(0, interval.getStart(), interval.getEnd() + 1);
            result.add(Triple.of(svInterval, interval, 0));
            final SimpleInterval simpleInterval10Padding = new SimpleInterval("seq1",
                    Math.max(1, interval.getStart() - 10), Math.min(100, interval.getEnd() + 10));
            result.add(Triple.of(svInterval, simpleInterval10Padding, 10));
        }
        result.add(Triple.of(new SVInterval(0, 1,2), new SimpleInterval("seq1", 1,1), 0));
        result.add(Triple.of(new SVInterval(0, 100, 101), new SimpleInterval("seq1", 100, 100), 0));
        result.add(Triple.of(new SVInterval(0, 100, 101), new SimpleInterval("seq1", 1, 100), 2000));

        result.add(Triple.of(new SVInterval(0, 1, 21), new SimpleInterval("seq1", 1, 40), 20));
        result.add(Triple.of(new SVInterval(2, 1, 10001), new SimpleInterval("seq3", 1, 10000), 0));
        result.add(Triple.of(new SVInterval(2, 1, 10001), new SimpleInterval("seq3", 1, 10000), 1));
        result.add(Triple.of(new SVInterval(2, 1, 10001), new SimpleInterval("seq3", 1, 10000), 1111));

        return result.stream()
                .map(t -> new Object[] {t.getLeft(), t.getMiddle(), t.getRight()})
                .toArray(Object[][]::new);
    }

    private static final SimpleInterval randomInterval(final Random rdn,
            final SAMSequenceDictionary dictionary,
            final String seqName, final int minStart, final int maxStart,
            final int minLength, final int maxLength) {
        final int seqIndex = seqName != null
                ? dictionary.getSequence(seqName).getSequenceIndex()
                : rdn.nextInt(dictionary.size());
        final int start = rdn.nextInt(maxStart - minStart) + minStart;
        final int length = rdn.nextInt(maxLength - minLength) + minLength;

        return new SimpleInterval(dictionary.getSequence(seqIndex).getSequenceName(),
                start, Math.min(length + start - 1, dictionary.getSequence(seqIndex).getSequenceLength()));
    }

}
