package org.broadinstitute.hellbender.tools.sv;

import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class SVDepthOnlyCallDefragmenterTest {

    private final static SVDepthOnlyCallDefragmenter defaultDefragmenter = new SVDepthOnlyCallDefragmenter(SVTestUtils.dict);

    private final static SVDepthOnlyCallDefragmenter singleSampleDefragmenter = new SVDepthOnlyCallDefragmenter(SVTestUtils.dict, SVDepthOnlyCallDefragmenter.DEFAULT_PADDING_FRACTION, 0, SVTestUtils.targetIntervals);

    @Test
    public void testFlattenCluster() {
        final SVCallRecordWithEvidence call1FlattenedDefault = defaultDefragmenter.flattenCluster(Collections.singletonList(SVTestUtils.call1));
        Assert.assertEquals(SVTestUtils.call1, call1FlattenedDefault);

        final SVCallRecordWithEvidence call1FlattenedSingleSample = singleSampleDefragmenter.flattenCluster(Collections.singletonList(SVTestUtils.call1));
        Assert.assertEquals(call1FlattenedSingleSample, call1FlattenedDefault);

        final SVCallRecordWithEvidence sameBoundsThreeSamples = singleSampleDefragmenter.flattenCluster(Arrays.asList(SVTestUtils.call1, SVTestUtils.sameBoundsSampleMismatch));
        Assert.assertEquals(sameBoundsThreeSamples.getStart(), SVTestUtils.call1.getStart());
        Assert.assertEquals(sameBoundsThreeSamples.getEnd(), SVTestUtils.call1.getEnd());
        Assert.assertTrue(sameBoundsThreeSamples.getGenotypes().get(0).sameGenotype(SVTestUtils.sample1));
        Assert.assertTrue(sameBoundsThreeSamples.getGenotypes().get(1).sameGenotype(SVTestUtils.sample2));
        Assert.assertTrue(sameBoundsThreeSamples.getGenotypes().get(2).sameGenotype(SVTestUtils.sample3));

        final SVCallRecordWithEvidence overlapping = singleSampleDefragmenter.flattenCluster(Arrays.asList(SVTestUtils.call1, SVTestUtils.call2));
        Assert.assertEquals(overlapping.getStart(), SVTestUtils.call1.getStart());
        Assert.assertEquals(overlapping.getEnd(), SVTestUtils.call2.getEnd());
    }

    @DataProvider
    public Object[][] clusterTogetherInputsDefault() {
        return new Object[][] {
                {SVTestUtils.call1, SVTestUtils.call1, true},
                {SVTestUtils.call1, SVTestUtils.call2, true},
                {SVTestUtils.call1, SVTestUtils.nonDepthOnly, false},
                {SVTestUtils.call1, SVTestUtils.sameBoundsSampleMismatch, false},
                {SVTestUtils.call1, SVTestUtils.nonDepthOnly, false}
        };
    }

    @DataProvider
    public Object[][] clusterTogetherInputsSingleSample() {
        return new Object[][] {
                {SVTestUtils.call1, SVTestUtils.call1, true},
                {SVTestUtils.call1, SVTestUtils.call2, true},  //overlapping, same samples
                {SVTestUtils.call1, SVTestUtils.nonDepthOnly, false},
                {SVTestUtils.call1, SVTestUtils.sameBoundsSampleMismatch, true},
                {SVTestUtils.call1, SVTestUtils.nonDepthOnly, false},
                {SVTestUtils.call1_CN1, SVTestUtils.call2_CN0, false}  //overlapping, but different copy number
        };
    }

    @Test(dataProvider = "clusterTogetherInputsDefault")
    public void testClusterTogetherDefault(final SVCallRecordWithEvidence call1, final SVCallRecordWithEvidence call2, final boolean expectedResult) {
        Assert.assertEquals(defaultDefragmenter.clusterTogether(call1, call2), expectedResult);
    }

    @Test(dataProvider = "clusterTogetherInputsSingleSample")
    public void testClusterTogetherSingleSample(final SVCallRecordWithEvidence call1, final SVCallRecordWithEvidence call2, final boolean expectedResult) {
        Assert.assertEquals(singleSampleDefragmenter.clusterTogether(call1, call2), expectedResult);
    }

    @Test
    public void testGetClusteringInterval() {
        Assert.assertTrue(defaultDefragmenter.getClusteringInterval(SVTestUtils.leftEdgeCall, null).getStart() > 0);
        Assert.assertTrue(singleSampleDefragmenter.getClusteringInterval(SVTestUtils.leftEdgeCall, null).getStart() > 0);
        Assert.assertEquals(defaultDefragmenter.getClusteringInterval(SVTestUtils.rightEdgeCall, null).getEnd(), SVTestUtils.chr1Length);
        Assert.assertTrue(singleSampleDefragmenter.getClusteringInterval(SVTestUtils.rightEdgeCall, null).getEnd() <= SVTestUtils.chr1Length);  //will be less than chr1length if target intervals are smaller than chr1


        final SimpleInterval littleCluster = new SimpleInterval("chr1", SVTestUtils.start, SVTestUtils.start + SVTestUtils.length -1);
        final SimpleInterval totalInterval = defaultDefragmenter.getClusteringInterval(SVTestUtils.call2, littleCluster);
        //interval describing cluster should already be padded
        Assert.assertEquals(totalInterval.getStart(), SVTestUtils.start);
        //padding is added to the input call
        Assert.assertEquals(totalInterval.getEnd(), SVTestUtils.call2.getEnd()+(int)Math.round(SVTestUtils.length*SVDepthOnlyCallDefragmenter.getDefaultPaddingFraction()));
    }

    @Test
    public void testAdd() {
        //single-sample merge case, ignoring sample sets
        final SVDepthOnlyCallDefragmenter temp1 = new SVDepthOnlyCallDefragmenter(SVTestUtils.dict, SVDepthOnlyCallDefragmenter.DEFAULT_PADDING_FRACTION, 0.8, SVTestUtils.targetIntervals);
        temp1.add(SVTestUtils.call1);
        //force new cluster by adding a non-overlapping event
        temp1.add(SVTestUtils.call3);
        final List<SVCallRecordWithEvidence> output1 = temp1.getOutput(); //flushes all clusters
        Assert.assertEquals(output1.size(), 2);
        Assert.assertEquals(SVTestUtils.call1, output1.get(0));
        Assert.assertEquals(SVTestUtils.call3, output1.get(1));

        final SVDepthOnlyCallDefragmenter temp2 = new SVDepthOnlyCallDefragmenter(SVTestUtils.dict, SVDepthOnlyCallDefragmenter.DEFAULT_PADDING_FRACTION, 0.8, SVTestUtils.targetIntervals);
        temp2.add(SVTestUtils.call1);
        temp2.add(SVTestUtils.call2);  //should overlap after padding
        //force new cluster by adding a call on another contig
        temp2.add(SVTestUtils.call4_chr10);
        final List<SVCallRecordWithEvidence> output2 = temp2.getOutput();
        Assert.assertEquals(output2.size(), 2);
        Assert.assertEquals(output2.get(0).getStart(), SVTestUtils.call1.getStart());
        Assert.assertEquals(output2.get(0).getEnd(), SVTestUtils.call2.getEnd());
        Assert.assertEquals(output2.get(1), SVTestUtils.call4_chr10);

        //cohort case, checking sample set overlap
        final SVDepthOnlyCallDefragmenter temp3 = new SVDepthOnlyCallDefragmenter(SVTestUtils.dict);
        temp3.add(SVTestUtils.call1);
        temp3.add(SVTestUtils.sameBoundsSampleMismatch);
        final List<SVCallRecordWithEvidence> output3 = temp3.getOutput();
        Assert.assertEquals(output3.size(), 2);
    }
}