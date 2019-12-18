package org.broadinstitute.hellbender.tools.sv;

import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.BeforeTest;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Collections;

public class SVDepthOnlyCallDefragmenterTest {

    final static SVDepthOnlyCallDefragmenter defaultDefragmenter = new SVDepthOnlyCallDefragmenter(SVTestUtils.dict);

    final static SVDepthOnlyCallDefragmenter singleSampleDefragmenter = new SVDepthOnlyCallDefragmenter(SVTestUtils.dict, 0.0);

    @BeforeTest
    public void initializeDefragmenters() {
        defaultDefragmenter.add(SVTestUtils.call1);
        singleSampleDefragmenter.add(SVTestUtils.call1);
    }

    @Test
    public void testFlattenCluster() {
        final SVCallRecordWithEvidence call1FlattenedDefault = defaultDefragmenter.flattenCluster(Collections.singletonList(SVTestUtils.call1));
        Assert.assertEquals(SVTestUtils.call1, call1FlattenedDefault);

        final SVCallRecordWithEvidence call1FlattenedSingleSample = singleSampleDefragmenter.flattenCluster(Collections.singletonList(SVTestUtils.call1));
        Assert.assertEquals(call1FlattenedSingleSample, call1FlattenedDefault);

        final SVCallRecordWithEvidence sameBoundsThreeSamples = singleSampleDefragmenter.flattenCluster(Arrays.asList(SVTestUtils.call1, SVTestUtils.sameBoundsSampleMismatch));
        Assert.assertEquals(sameBoundsThreeSamples.getStart(), SVTestUtils.call1.getStart());
        Assert.assertEquals(sameBoundsThreeSamples.getEnd(), SVTestUtils.call1.getEnd());
        Assert.assertTrue(sameBoundsThreeSamples.getGenotypes().containsAll(SVTestUtils.threeGenotypes));

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
                {SVTestUtils.call1, SVTestUtils.call2, true},
                {SVTestUtils.call1, SVTestUtils.nonDepthOnly, false},
                {SVTestUtils.call1, SVTestUtils.sameBoundsSampleMismatch, true},
                {SVTestUtils.call1, SVTestUtils.nonDepthOnly, false}
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
        Assert.assertTrue(defaultDefragmenter.getClusteringInterval(SVTestUtils.rightEdgeCall, null).getEnd() == SVTestUtils.chr1Length);
        Assert.assertTrue(singleSampleDefragmenter.getClusteringInterval(SVTestUtils.rightEdgeCall, null).getEnd() == SVTestUtils.chr1Length);


        final SimpleInterval littleCluster = new SimpleInterval("chr1", SVTestUtils.start, SVTestUtils.start + SVTestUtils.length -1);
        final SimpleInterval totalInterval = defaultDefragmenter.getClusteringInterval(SVTestUtils.call2, littleCluster);
        //interval describing cluster should already be padded
        Assert.assertEquals(totalInterval.getStart(), SVTestUtils.start);
        //padding is added to the input call
        Assert.assertEquals(totalInterval.getEnd(), SVTestUtils.call2.getEnd()+(int)Math.round(SVTestUtils.length*defaultDefragmenter.getPaddingFraction()));
    }
}