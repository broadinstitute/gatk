package org.broadinstitute.hellbender.tools.sv.cluster;

import com.google.common.collect.Lists;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
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
    private static final SVClusterEngine<SVCallRecord> defaultDefragmenter = SVClusterEngineFactory.createCNVDefragmenter(SVTestUtils.hg38Dict, CanonicalSVCollapser.AltAlleleSummaryStrategy.COMMON_SUBTYPE, SVTestUtils.hg38Reference, paddingFraction, sampleOverlap, CanonicalSVLinkage.DEFAULT_DEPTH_ONLY_PARAMS);
    private static final SVClusterEngine<SVCallRecord> binnedDefragmenter = SVClusterEngineFactory.createBinnedCNVDefragmenter(SVTestUtils.hg38Dict, CanonicalSVCollapser.AltAlleleSummaryStrategy.COMMON_SUBTYPE, SVTestUtils.hg38Reference, paddingFraction, 0, SVTestUtils.targetIntervals, CanonicalSVLinkage.DEFAULT_DEPTH_ONLY_PARAMS);

    @Test
    public void testCollapser() {
        final SVCallRecord call1FlattenedDefault = defaultDefragmenter.getCollapser().collapse(Collections.singletonList(SVTestUtils.call1));
        SVTestUtils.assertEqualsExceptMembershipAndGT(SVTestUtils.call1, call1FlattenedDefault);

        final SVCallRecord call1FlattenedSingleSample = binnedDefragmenter.getCollapser().collapse(Collections.singletonList(SVTestUtils.call1));
        SVTestUtils.assertEqualsExceptMembershipAndGT(call1FlattenedSingleSample, call1FlattenedDefault);

        final SVCallRecord sameBoundsThreeSamples = binnedDefragmenter.getCollapser().collapse(Arrays.asList(SVTestUtils.call1, SVTestUtils.sameBoundsSampleMismatch));
        Assert.assertEquals(sameBoundsThreeSamples.getPositionA(), SVTestUtils.call1.getPositionA());
        Assert.assertEquals(sameBoundsThreeSamples.getPositionB(), SVTestUtils.call1.getPositionB());

        final Genotype testGenotype1 = sameBoundsThreeSamples.getGenotypes().get(SVTestUtils.sample1.make().getSampleName());
        final Genotype expectedGenotype1 = SVTestUtils.sample1.alleles(Lists.newArrayList(Allele.REF_T, Allele.SV_SIMPLE_DEL)).make();
        Assert.assertEquals(testGenotype1.getAlleles(), expectedGenotype1.getAlleles());
        Assert.assertEquals(testGenotype1.getExtendedAttributes(), expectedGenotype1.getExtendedAttributes());

        final Genotype testGenotype2 = sameBoundsThreeSamples.getGenotypes().get(SVTestUtils.sample2.make().getSampleName());
        final Genotype expectedGenotype2 = SVTestUtils.sample2.alleles(Collections.nCopies(2, Allele.NO_CALL)).make();
        Assert.assertEquals(testGenotype2.getAlleles(), expectedGenotype2.getAlleles());
        Assert.assertEquals(testGenotype2.getExtendedAttributes(), expectedGenotype2.getExtendedAttributes());

        final Genotype testGenotype3 = sameBoundsThreeSamples.getGenotypes().get(SVTestUtils.sample3.make().getSampleName());
        final Genotype expectedGenotype3 = SVTestUtils.sample3.alleles(Lists.newArrayList(Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DEL)).make();
        Assert.assertEquals(testGenotype3.getAlleles(), expectedGenotype3.getAlleles());
        Assert.assertEquals(testGenotype3.getExtendedAttributes(), expectedGenotype3.getExtendedAttributes());

        final SVCallRecord overlapping = binnedDefragmenter.getCollapser().collapse(Arrays.asList(SVTestUtils.call1, SVTestUtils.call2));
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
        Assert.assertEquals(defaultDefragmenter.getLinkage().areClusterable(call1, call2), expectedResult, name);
    }

    @Test(dataProvider = "clusterTogetherInputsSingleSample")
    public void testClusterTogetherSingleSample(final SVCallRecord call1, final SVCallRecord call2,
                                                final boolean expectedResult, final String name) {
        Assert.assertEquals(binnedDefragmenter.getLinkage().areClusterable(call1, call2), expectedResult, name);
    }

    @Test
    public void testGetMaxClusterableStartingPosition() {
        Assert.assertEquals(defaultDefragmenter.getLinkage().getMaxClusterableStartingPosition(SVTestUtils.rightEdgeCall), SVTestUtils.chr1Length);
        Assert.assertTrue(binnedDefragmenter.getLinkage().getMaxClusterableStartingPosition(SVTestUtils.rightEdgeCall) == SVTestUtils.chr1Length);  //will be less than chr1length if target intervals are smaller than chr1
    }

    @Test
    public void testAdd() {
        //single-sample merge case, ignoring sample sets
        final SVClusterEngine<SVCallRecord> temp1 = SVClusterEngineFactory.createBinnedCNVDefragmenter(SVTestUtils.hg38Dict, CanonicalSVCollapser.AltAlleleSummaryStrategy.COMMON_SUBTYPE, SVTestUtils.hg38Reference, paddingFraction, 0.8, SVTestUtils.targetIntervals, CanonicalSVLinkage.DEFAULT_DEPTH_ONLY_PARAMS);
        temp1.add(SVTestUtils.call1);
        //force new cluster by adding a non-overlapping event
        temp1.add(SVTestUtils.call3);
        final List<SVCallRecord> output1 = temp1.getOutput(); //flushes all clusters
        Assert.assertEquals(output1.size(), 2);
        SVTestUtils.assertEqualsExceptMembershipAndGT(SVTestUtils.call1, output1.get(0));
        SVTestUtils.assertEqualsExceptMembershipAndGT(SVTestUtils.call3, output1.get(1));

        final SVClusterEngine<SVCallRecord> temp2 = SVClusterEngineFactory.createBinnedCNVDefragmenter(SVTestUtils.hg38Dict, CanonicalSVCollapser.AltAlleleSummaryStrategy.COMMON_SUBTYPE, SVTestUtils.hg38Reference, paddingFraction, 0.8, SVTestUtils.targetIntervals, CanonicalSVLinkage.DEFAULT_DEPTH_ONLY_PARAMS);
        temp2.add(SVTestUtils.call1);
        temp2.add(SVTestUtils.call2);  //should overlap after padding
        //force new cluster by adding a call on another contig
        temp2.add(SVTestUtils.call4_chr10);
        final List<SVCallRecord> output2 = temp2.getOutput();
        Assert.assertEquals(output2.size(), 2);
        Assert.assertEquals(output2.get(0).getPositionA(), SVTestUtils.call1.getPositionA());
        Assert.assertEquals(output2.get(0).getPositionB(), SVTestUtils.call2.getPositionB());
        SVTestUtils.assertEqualsExceptMembershipAndGT(output2.get(1), SVTestUtils.call4_chr10);

        //cohort case, checking sample set overlap
        final SVClusterEngine<SVCallRecord> temp3 = SVClusterEngineFactory.createCNVDefragmenter(SVTestUtils.hg38Dict, CanonicalSVCollapser.AltAlleleSummaryStrategy.COMMON_SUBTYPE, SVTestUtils.hg38Reference, CNVLinkage.DEFAULT_PADDING_FRACTION, CNVLinkage.DEFAULT_SAMPLE_OVERLAP, CanonicalSVLinkage.DEFAULT_DEPTH_ONLY_PARAMS);
        temp3.add(SVTestUtils.call1);
        temp3.add(SVTestUtils.sameBoundsSampleMismatch);
        final List<SVCallRecord> output3 = temp3.getOutput();
        Assert.assertEquals(output3.size(), 2);
    }
}