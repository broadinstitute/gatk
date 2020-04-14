package org.broadinstitute.hellbender.tools.sv;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.StructuralVariantType;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.BeforeTest;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import sun.nio.ch.SelectorImpl;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class SVDepthOnlyCallDefragmenterTest {

    final static int chr1Length = 249250621;

    final static SAMSequenceDictionary dict = new SAMSequenceDictionary(
            Arrays.asList(new SAMSequenceRecord("chr1", chr1Length),
                    new SAMSequenceRecord("chr10", 135534747),
                    new SAMSequenceRecord("chrX", 155270560)));

    final static SVDepthOnlyCallDefragmenter defaultDefragmenter = new SVDepthOnlyCallDefragmenter(dict);

    final static SVDepthOnlyCallDefragmenter singleSampleDefragmenter = new SVDepthOnlyCallDefragmenter(dict, 0.0);

    final static int start = 10001;
    final static int length = 10000;
    final static int start2 = (start + length -1) + (int)Math.round(length * defaultDefragmenter.getPaddingFraction());

    final static Genotype sample1 = GenotypeBuilder.create("sample1", Collections.singletonList(Allele.create("<"+GATKSVVCFConstants.SYMB_ALT_ALLELE_DEL+">", false)));

    final static Genotype sample2 = GenotypeBuilder.create("sample2", Collections.singletonList(Allele.create("<"+GATKSVVCFConstants.SYMB_ALT_ALLELE_DUP+">", false)));

    final static Genotype sample3 = GenotypeBuilder.create("sample3", Collections.singletonList(Allele.create("<"+GATKSVVCFConstants.SYMB_ALT_ALLELE_DEL+">", false)));

    final static List<Genotype> threeGenotypes = Arrays.asList(sample1, sample2, sample3);

    final static SVCallRecordWithEvidence call1 = new SVCallRecordWithEvidence("chr1", start, true,
            "chr1", 20000, true,
            StructuralVariantType.CNV, 10000,
            Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
            Arrays.asList(sample1, sample2),
            Collections.emptyList(), Collections.emptyList(), Collections.emptyList());

    final static SVCallRecordWithEvidence call2 = new SVCallRecordWithEvidence("chr1", start2, true,
            "chr1", start2 + length -1, true,
            StructuralVariantType.CNV, length,
            Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
            Arrays.asList(sample1, sample2),
            Collections.emptyList(), Collections.emptyList(), Collections.emptyList());

    final static SVCallRecordWithEvidence sameBoundsSampleMismatch = new SVCallRecordWithEvidence("chr1", start, true,
            "chr1", start + length -1, true,
            StructuralVariantType.CNV, length,
            Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
            Arrays.asList(sample1, sample3),
            Collections.emptyList(), Collections.emptyList(), Collections.emptyList());

    final static SVCallRecordWithEvidence nonDepthOnly = new SVCallRecordWithEvidence("chr1", start, true,
            "chr1", 20000, true,
            StructuralVariantType.CNV, length,
            Arrays.asList(GATKSVVCFConstants.DEPTH_ALGORITHM, "PE"),
            Arrays.asList(sample1, sample2),
            Collections.emptyList(), Collections.emptyList(), Collections.emptyList());

    final static SVCallRecordWithEvidence leftEdgeCall = new SVCallRecordWithEvidence("chr1", 1, true,
            "chr1", length, true,
            StructuralVariantType.CNV, length,
            Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
            Arrays.asList(sample1, sample2),
            Collections.emptyList(), Collections.emptyList(), Collections.emptyList());

    final static SVCallRecordWithEvidence rightEdgeCall = new SVCallRecordWithEvidence("chr1", chr1Length - 99, true,
            "chr1", chr1Length, true,
            StructuralVariantType.CNV, 100,
            Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
            Arrays.asList(sample1, sample2),
            Collections.emptyList(), Collections.emptyList(), Collections.emptyList());

    @BeforeTest
    public void initializeDefragmenters() {
        defaultDefragmenter.add(call1);
        singleSampleDefragmenter.add(call1);
    }

    @Test
    public void testFlattenCluster() {
        final SVCallRecordWithEvidence call1FlattenedDefault = defaultDefragmenter.flattenCluster(Collections.singletonList(call1));
        Assert.assertEquals(call1, call1FlattenedDefault);

        final SVCallRecordWithEvidence call1FlattenedSingleSample = singleSampleDefragmenter.flattenCluster(Collections.singletonList(call1));
        Assert.assertEquals(call1FlattenedSingleSample, call1FlattenedDefault);

        final SVCallRecordWithEvidence sameBoundsThreeSamples = singleSampleDefragmenter.flattenCluster(Arrays.asList(call1, sameBoundsSampleMismatch));
        Assert.assertEquals(sameBoundsThreeSamples.getStart(), call1.getStart());
        Assert.assertEquals(sameBoundsThreeSamples.getEnd(), call1.getEnd());
        Assert.assertTrue(sameBoundsThreeSamples.getGenotypes().containsAll(threeGenotypes));

        final SVCallRecordWithEvidence overlapping = singleSampleDefragmenter.flattenCluster(Arrays.asList(call1, call2));
        Assert.assertEquals(overlapping.getStart(), call1.getStart());
        Assert.assertEquals(overlapping.getEnd(), call2.getEnd());
    }

    @DataProvider
    public Object[][] clusterTogetherInputsDefault() {
        return new Object[][] {
                {call1, call1, true},
                {call1, call2, true},
                {call1, nonDepthOnly, false},
                {call1, sameBoundsSampleMismatch, false},
                {call1, nonDepthOnly, false}
        };
    }

    @DataProvider
    public Object[][] clusterTogetherInputsSingleSample() {
        return new Object[][] {
                {call1, call1, true},
                {call1, call2, true},
                {call1, nonDepthOnly, false},
                {call1, sameBoundsSampleMismatch, true},
                {call1, nonDepthOnly, false}
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
        Assert.assertTrue(defaultDefragmenter.getClusteringInterval(leftEdgeCall, null).getStart() > 0);
        Assert.assertTrue(singleSampleDefragmenter.getClusteringInterval(leftEdgeCall, null).getStart() > 0);
        Assert.assertTrue(singleSampleDefragmenter.getClusteringInterval(rightEdgeCall, null).getEnd() == chr1Length);
        Assert.assertTrue(singleSampleDefragmenter.getClusteringInterval(rightEdgeCall, null).getEnd() == chr1Length);


        final SimpleInterval littleCluster = new SimpleInterval("chr1", start, start + length -1);
        final SimpleInterval totalInterval = defaultDefragmenter.getClusteringInterval(call2, littleCluster);
        //interval describing cluster should already be padded
        Assert.assertEquals(totalInterval.getStart(), start);
        //padding is added to the input call
        Assert.assertEquals(totalInterval.getEnd(), call2.getEnd()+(int)Math.round(length*defaultDefragmenter.getPaddingFraction()));
    }
}