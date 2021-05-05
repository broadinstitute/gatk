package org.broadinstitute.hellbender.tools.sv.cluster;

import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SVTestUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class BinnedCNVDefragmenterTest {

    private static final double paddingFraction = 0.5;
    private static final double sampleOverlap = 0.9;
    private static final CNVDefragmenter defaultDefragmenter = new CNVDefragmenter(SVTestUtils.dict, paddingFraction, sampleOverlap);
    private static final BinnedCNVDefragmenter binnedDefragmenter = new BinnedCNVDefragmenter(SVTestUtils.dict, paddingFraction, 0, SVTestUtils.targetIntervals);

    @Test
    public void testCollapser() {
        final SVCallRecord call1FlattenedDefault = defaultDefragmenter.getCollapser().apply(Collections.singletonList(SVTestUtils.call1));
        SVTestUtils.assertEqualsExceptMembership(SVTestUtils.call1, call1FlattenedDefault);

        final SVCallRecord call1FlattenedSingleSample = binnedDefragmenter.getCollapser().apply(Collections.singletonList(SVTestUtils.call1));
        SVTestUtils.assertEqualsExceptMembership(call1FlattenedSingleSample, call1FlattenedDefault);

        final SVCallRecord sameBoundsThreeSamples = binnedDefragmenter.getCollapser().apply(Arrays.asList(SVTestUtils.call1, SVTestUtils.sameBoundsSampleMismatch));
        Assert.assertEquals(sameBoundsThreeSamples.getPositionA(), SVTestUtils.call1.getPositionA());
        Assert.assertEquals(sameBoundsThreeSamples.getPositionB(), SVTestUtils.call1.getPositionB());
        Assert.assertTrue(sameBoundsThreeSamples.getGenotypes().get(SVTestUtils.sample1.getSampleName()).sameGenotype(SVTestUtils.sample1));
        Assert.assertTrue(sameBoundsThreeSamples.getGenotypes().get(SVTestUtils.sample2.getSampleName()).sameGenotype(SVTestUtils.sample2));
        Assert.assertTrue(sameBoundsThreeSamples.getGenotypes().get(SVTestUtils.sample3.getSampleName()).sameGenotype(SVTestUtils.sample3));

        final SVCallRecord overlapping = binnedDefragmenter.getCollapser().apply(Arrays.asList(SVTestUtils.call1, SVTestUtils.call2));
        Assert.assertEquals(overlapping.getPositionA(), SVTestUtils.call1.getPositionA());
        Assert.assertEquals(overlapping.getPositionB(), SVTestUtils.call2.getPositionB());
    }

    @DataProvider
    public Object[][] clusterTogetherInputsDefault() {
        return new Object[][] {
                {SVTestUtils.call1, SVTestUtils.call1, true, "call1 call1"},
                {SVTestUtils.call1, SVTestUtils.call2, true, "call1 call2"},
                {SVTestUtils.call1, SVTestUtils.nonDepthOnly, false, "call1 nonDepthOnly"},
                {SVTestUtils.call1, SVTestUtils.sameBoundsSampleMismatch, false, "call1 sameBoundsSampleMismatch"}
        };
    }

    @DataProvider
    public Object[][] clusterTogetherInputsSingleSample() {
        return new Object[][] {
                {SVTestUtils.call1, SVTestUtils.call1, true, "call1 call1"},
                {SVTestUtils.call1, SVTestUtils.call2, true, "call1 call2"},  //overlapping, same samples
                {SVTestUtils.call1, SVTestUtils.nonDepthOnly, false, "call1 nonDepthOnly"},
                {SVTestUtils.call1, SVTestUtils.sameBoundsSampleMismatch, true, "call1 sameBoundsSampleMismatch"},
                {SVTestUtils.call1_CN1, SVTestUtils.call2_CN0, false, "call1_CN1 call2_CN0"}  //overlapping, but different copy number
        };
    }

    @Test(dataProvider = "clusterTogetherInputsDefault")
    public void testClusterTogetherDefault(final SVCallRecord call1, final SVCallRecord call2,
                                           final boolean expectedResult, final String name) {
        Assert.assertEquals(defaultDefragmenter.clusterTogether(call1, call2), expectedResult, name);
    }

    @Test(dataProvider = "clusterTogetherInputsSingleSample")
    public void testClusterTogetherSingleSample(final SVCallRecord call1, final SVCallRecord call2,
                                                final boolean expectedResult, final String name) {
        Assert.assertEquals(binnedDefragmenter.clusterTogether(call1, call2), expectedResult, name);
    }

    @Test
    public void testGetMaxClusterableStartingPosition() {
        Assert.assertEquals(defaultDefragmenter.getMaxClusterableStartingPosition(SVTestUtils.rightEdgeCall), SVTestUtils.chr1Length);
        Assert.assertTrue(binnedDefragmenter.getMaxClusterableStartingPosition(SVTestUtils.rightEdgeCall) == SVTestUtils.chr1Length);  //will be less than chr1length if target intervals are smaller than chr1
    }

    @Test
    public void testAdd() {
        //single-sample merge case, ignoring sample sets
        final CNVDefragmenter temp1 = new BinnedCNVDefragmenter(SVTestUtils.dict, paddingFraction, 0.8, SVTestUtils.targetIntervals);
        temp1.add(SVTestUtils.call1);
        //force new cluster by adding a non-overlapping event
        temp1.add(SVTestUtils.call3);
        final List<SVCallRecord> output1 = temp1.getOutput(); //flushes all clusters
        Assert.assertEquals(output1.size(), 2);
        SVTestUtils.assertEqualsExceptMembership(SVTestUtils.call1, output1.get(0));
        SVTestUtils.assertEqualsExceptMembership(SVTestUtils.call3, output1.get(1));

        final CNVDefragmenter temp2 = new BinnedCNVDefragmenter(SVTestUtils.dict, paddingFraction, 0.8, SVTestUtils.targetIntervals);
        temp2.add(SVTestUtils.call1);
        temp2.add(SVTestUtils.call2);  //should overlap after padding
        //force new cluster by adding a call on another contig
        temp2.add(SVTestUtils.call4_chr10);
        final List<SVCallRecord> output2 = temp2.getOutput();
        Assert.assertEquals(output2.size(), 2);
        Assert.assertEquals(output2.get(0).getPositionA(), SVTestUtils.call1.getPositionA());
        Assert.assertEquals(output2.get(0).getPositionB(), SVTestUtils.call2.getPositionB());
        SVTestUtils.assertEqualsExceptMembership(output2.get(1), SVTestUtils.call4_chr10);

        //cohort case, checking sample set overlap
        final CNVDefragmenter temp3 = new CNVDefragmenter(SVTestUtils.dict, CNVDefragmenter.DEFAULT_PADDING_FRACTION, CNVDefragmenter.DEFAULT_SAMPLE_OVERLAP);
        temp3.add(SVTestUtils.call1);
        temp3.add(SVTestUtils.sameBoundsSampleMismatch);
        final List<SVCallRecord> output3 = temp3.getOutput();
        Assert.assertEquals(output3.size(), 2);
    }
}