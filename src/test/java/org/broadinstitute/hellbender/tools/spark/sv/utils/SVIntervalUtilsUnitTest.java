package org.broadinstitute.hellbender.tools.spark.sv.utils;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.EvidenceTargetLink;
import org.broadinstitute.hellbender.utils.GenomeLoc;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

public class SVIntervalUtilsUnitTest extends BaseTest {

    private static EvidenceTargetLink emptyLink(final int contig) {
        final StrandedInterval source = new StrandedInterval(new SVInterval(contig, 100, 200), true);
        final StrandedInterval target = new StrandedInterval(new SVInterval(contig, 1100, 1200), false);
        return new EvidenceTargetLink(source, target, 0, 0, Collections.emptySet(), Collections.emptySet());
    }

    @Test(groups = "sv")
    public void getOuterIntrachromosomalLinkIntervalTest() {
        final StrandedInterval source = new StrandedInterval(new SVInterval(0, 100, 200), true);
        final StrandedInterval target = new StrandedInterval(new SVInterval(0, 1100, 1200), false);
        final Set<String> readPairTemplates = new HashSet<>();
        readPairTemplates.add("READ_1");
        readPairTemplates.add("READ_2");
        final EvidenceTargetLink link = new EvidenceTargetLink(source, target, 0, readPairTemplates.size(), readPairTemplates, Collections.emptySet());
        final SVInterval linkOuterInterval = SVIntervalUtils.getOuterIntrachromosomalLinkInterval(link);
        Assert.assertEquals(0, linkOuterInterval.getContig());
        Assert.assertEquals(100, linkOuterInterval.getStart());
        Assert.assertEquals(1200, linkOuterInterval.getEnd());
    }

    @Test(groups = "sv")
    public void getPaddedIntervalTest() {
        final SVInterval interval = new SVInterval(0, 100, 200);
        final int padding = 10;
        final SAMSequenceDictionary dictionary = new SAMSequenceDictionary();
        dictionary.addSequence(new SAMSequenceRecord("seq0", 1000));
        final SVInterval paddedInterval = SVIntervalUtils.getPaddedInterval(interval, padding, dictionary);
        Assert.assertEquals(paddedInterval.getContig(), 0);
        Assert.assertEquals(paddedInterval.getStart(), 90);
        Assert.assertEquals(paddedInterval.getEnd(), 210);

        final SVInterval interval2 = new SVInterval(0, 5, 995);
        final SVInterval paddedInterval2 = SVIntervalUtils.getPaddedInterval(interval2, padding, dictionary);
        Assert.assertEquals(paddedInterval2.getContig(), 0);
        Assert.assertEquals(paddedInterval2.getStart(), 0);
        Assert.assertEquals(paddedInterval2.getEnd(), 999);
    }

    @Test(groups = "sv")
    public void intervalTreeSearchTest() {
        final SVIntervalTree<EvidenceTargetLink> tree = new SVIntervalTree<>();
        final EvidenceTargetLink link0 = emptyLink(0);
        final EvidenceTargetLink link1 = emptyLink(1);
        tree.put(new SVInterval(0, 100, 200), link0);
        tree.put(new SVInterval(0, 300, 400), link1);
        final SVInterval interval = new SVInterval(1, 0, 500);

        final Collection<EvidenceTargetLink> emptyOverlapResult = SVIntervalUtils.getOverlappingLinksOnInterval(interval, tree);
        Assert.assertTrue(emptyOverlapResult.isEmpty());
        Assert.assertFalse(SVIntervalUtils.hasReciprocalOverlapInTree(interval, tree, 0));

        final SVInterval interval2 = new SVInterval(0, 150, 250);
        final Collection<EvidenceTargetLink> overlapResult = SVIntervalUtils.getOverlappingLinksOnInterval(interval2, tree);
        final Iterator<EvidenceTargetLink> resultIter = overlapResult.iterator();
        Assert.assertTrue(resultIter.hasNext());
        Assert.assertEquals(resultIter.next(), link0);
        Assert.assertEquals(overlapResult.size(), 1);

        Assert.assertFalse(SVIntervalUtils.hasReciprocalOverlapInTree(interval2, tree, 0.6));
        Assert.assertTrue(SVIntervalUtils.hasReciprocalOverlapInTree(interval2, tree, 0.5));
    }

    @DataProvider(name = "containsIntervalData")
    public Object[][] getContainsIntervalData() {
        return new Object[][]{
                {0, 10, 11, 12, false},
                {0, 10, 10, 12, false},
                {0, 10, 0, 10, true},
                {0, 10, 0, 9, true},
                {0, 10, 1, 10, true},
                {0, 10, 4, 8, true},
                {0, 10, -1, 0, false},
                {0, 10, -1, 1, false},
                {0, 10, -1, 10, false},
                {0, 10, -1, 11, false},
                {0, 10, 11, 15, false},
                {0, 10, 0, 0, true},
                {0, 10, 10, 10, true},
                {0, 10, 0, 0, true},
                {0, 10, 5, 5, true}
        };
    }

    @Test(groups = "sv",
            dataProvider = "containsIntervalData")
    public void containsIntervalTest(final int startA, final int endA, final int startB, final int endB, final boolean expectedResult) {
        final SVInterval intervalA = new SVInterval(0, startA, endA);
        final SVInterval intervalB1 = new SVInterval(0, startB, endB);
        Assert.assertEquals(SVIntervalUtils.containsInterval(intervalA, intervalB1), expectedResult);
        final SVInterval intervalB2 = new SVInterval(1, startB, endB);
        Assert.assertEquals(SVIntervalUtils.containsInterval(intervalA, intervalB2), false);
    }

    @DataProvider(name = "reciprocalOverlapTest")
    public Object[][] getReciprocalOverlapData() {
        return new Object[][]{
                {10, 20, 0, 9, 0., 1.},
                {10, 20, 20, 29, 0., 1.},
                {10, 20, 5, 15, 5. / 10., 6. / 10.},
                {10, 20, 15, 25, 5. / 10., 6. / 10.},
                {10, 20, 19, 20, 1. / 10., 2. / 10.},
                {10, 20, 19, 21, 1. / 10., 2. / 10.},
                {10, 20, 18, 21, 2. / 10., 3. / 10.},
                {10, 20, 18, 30, 2. / 12., 3. / 12.}
        };
    }

    @Test(groups = "sv",
            dataProvider = "reciprocalOverlapTest")
    public void reciprocalOverlapTest(final int startA, final int endA, final int startB, final int endB, final double reciprocalOverlapEquals, final double reciprocalOverlapFalse) {
        final SVInterval intervalA = new SVInterval(0, startA, endA);
        final SVInterval intervalB1 = new SVInterval(0, startB, endB);
        Assert.assertEquals(SVIntervalUtils.reciprocalOverlap(intervalA, intervalB1), reciprocalOverlapEquals);
        Assert.assertTrue(SVIntervalUtils.hasReciprocalOverlap(intervalA, intervalB1, reciprocalOverlapEquals));
        Assert.assertTrue(SVIntervalUtils.hasReciprocalOverlap(intervalA, intervalB1, reciprocalOverlapEquals / 2));
        Assert.assertFalse(SVIntervalUtils.hasReciprocalOverlap(intervalA, intervalB1, reciprocalOverlapFalse));
        final SVInterval intervalB2 = new SVInterval(1, startB, endB);
        Assert.assertEquals(SVIntervalUtils.reciprocalOverlap(intervalA, intervalB2), 0.);
    }

    @Test(groups = "sv")
    public void convertIntervalsTest() {
        final SVInterval svInterval = new SVInterval(0, 100, 200);
        final SAMSequenceDictionary dictionary = new SAMSequenceDictionary();
        dictionary.addSequence(new SAMSequenceRecord("seq0", 1000));
        final SimpleInterval simpleInterval = SVIntervalUtils.convertToSimpleInterval(svInterval, dictionary);
        Assert.assertEquals(simpleInterval, new SimpleInterval("seq0", 100, 200));
        final SVInterval convertedSVInterval = SVIntervalUtils.convertToSVInterval(simpleInterval, dictionary);
        Assert.assertEquals(convertedSVInterval, svInterval);
        final GenomeLoc genomeLoc = SVIntervalUtils.convertToGenomeLoc(svInterval, dictionary);
        Assert.assertEquals(genomeLoc, new GenomeLoc("seq0", 0, 100, 200));
    }

}