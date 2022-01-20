package org.broadinstitute.hellbender.tools.gvs.common;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.OverlapDetector;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.TreeSet;


public class WeightedSplitIntervalsTest {
    private SAMSequenceDictionary header =
            new SAMSequenceDictionary(Collections.singletonList(
                    new SAMSequenceRecord("chr1", 200000000)));

    @Test
    public void testTwoWeights() {
        Interval i = new Interval("chr1", 1000, 1999);
        List<WeightedInterval> weights = new ArrayList<>();
        weights.add(new WeightedInterval("chr1", 1000, 1499, 500));
        weights.add(new WeightedInterval("chr1", 1500, 1999, 1000));
        OverlapDetector<WeightedInterval> od = OverlapDetector.create(weights);

        TreeSet<WeightedInterval> out = new TreeSet<>(WeightedSplitIntervals.applyWeightsToInterval(header, i, od, 0));

        assertWeight(out.pollFirst(), "chr1", 1000, 1499, 500);
        assertWeight(out.pollFirst(), "chr1", 1500, 1999, 1000);

    }

    @Test
    public void testHalfWeight() {
        Interval i = new Interval("chr1", 1000, 1999);
        List<WeightedInterval> weights = new ArrayList<>();
        weights.add(new WeightedInterval("chr1", 1500, 2499, 4000));
        OverlapDetector<WeightedInterval> od = OverlapDetector.create(weights);

        TreeSet<WeightedInterval> out = new TreeSet<>(WeightedSplitIntervals.applyWeightsToInterval(header, i, od, 0));

        assertWeight(out.pollFirst(), "chr1", 1000, 1499, 0);
        assertWeight(out.pollFirst(), "chr1", 1500, 1999, 2000);
    }

    @Test
    public void testComplexWeight() {
        // complex test with weight halfway over beginning, no weight over middle, different weight over more bases, finishing with no weight at end
        Interval i = new Interval("chr1", 1000, 1999);

        List<WeightedInterval> weights = new ArrayList<>();
        weights.add(new WeightedInterval("chr1", 500, 1499, 4000)); // 4 per base
        weights.add(new WeightedInterval("chr1", 1600, 1699, 500)); // 5 per base
        OverlapDetector<WeightedInterval> od = OverlapDetector.create(weights);

        TreeSet<WeightedInterval> out = new TreeSet<>(WeightedSplitIntervals.applyWeightsToInterval(header, i, od, 0));

        assertWeight(out.pollFirst(), "chr1", 1000, 1499, 2000);
        assertWeight(out.pollFirst(), "chr1", 1500, 1599, 0);
        assertWeight(out.pollFirst(), "chr1", 1600, 1699, 500);
        assertWeight(out.pollFirst(), "chr1", 1700, 1999, 0);
    }

    private void assertWeight(WeightedInterval w, String contig, int start, int end, long weight) {
        Assert.assertEquals(w.getContig(), contig);
        Assert.assertEquals(w.getStart(), start);
        Assert.assertEquals(w.getEnd(), end);
        Assert.assertEquals(w.getWeight(), weight);
    }
}


