package org.broadinstitute.hellbender.tools.sv.cluster;

import com.google.common.collect.Lists;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.StructuralVariantType;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SVTestUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Collections;

public class CNVDefragmenterTest {

    private final SAMSequenceDictionary dictionary = SVTestUtils.dict;
    private final CNVDefragmenter defragmenter = new CNVDefragmenter(dictionary, 0.5, 0.6);

    @Test
    public void testClusterTogether() {
        final SVCallRecord deletion = new SVCallRecord("test_del", "chr1", 1000, true, "chr1", 1999, false, StructuralVariantType.DEL,
                1000, Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL),
                Collections.emptyList(), Collections.emptyMap());
        final SVCallRecord duplication = new SVCallRecord("test_dup", "chr1", 1000, false, "chr1", 1999, false, StructuralVariantType.DUP,
                1000, Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DUP),
                Collections.emptyList(), Collections.emptyMap());
        Assert.assertFalse(defragmenter.clusterTogether(deletion, duplication), "Different sv types should not cluster");

        final SVCallRecord duplicationNonDepthOnly = new SVCallRecord("test_dup", "chr1", 1000, false, "chr1", 1999, false, StructuralVariantType.DUP,
                1000, Lists.newArrayList(GATKSVVCFConstants.DEPTH_ALGORITHM, "pesr"),
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DUP),
                Collections.emptyList(), Collections.emptyMap());
        Assert.assertFalse(defragmenter.clusterTogether(duplication, duplicationNonDepthOnly), "Clustered records must be depth-only");

        final SVCallRecord cnv = new SVCallRecord("test_cnv", "chr1", 1000, false, "chr1", 1999, false, StructuralVariantType.CNV,
                1000, Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DUP),
                Collections.emptyList(), Collections.emptyMap());
        Assert.assertFalse(defragmenter.clusterTogether(deletion, cnv), "Different sv types should not cluster");

        final SVCallRecord insertion = new SVCallRecord("test_ins", "chr1", 1000, true, "chr1", 1001, false, StructuralVariantType.INS,
                1000, Collections.singletonList("pesr"),
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_INS),
                Collections.emptyList(), Collections.emptyMap());
        Assert.assertFalse(defragmenter.clusterTogether(insertion, insertion), "Only CNVs should be valid");

        final SVCallRecord deletion2 = new SVCallRecord("test_del2", "chr1", 1000, true, "chr1", 1999, false, StructuralVariantType.DEL,
                1000, Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL),
                Collections.emptyList(), Collections.emptyMap());
        Assert.assertTrue(defragmenter.clusterTogether(deletion, deletion2), "Valid identical records should cluster");

        final SVCallRecord deletion3 = new SVCallRecord("test_del3", "chr1", 2999, true, "chr1", 3998, false, StructuralVariantType.DEL,
                1000, Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL),
                Collections.emptyList(), Collections.emptyMap());
        Assert.assertTrue(defragmenter.clusterTogether(deletion, deletion3), "Should cluster due to overlap");

        final SVCallRecord deletion4 = new SVCallRecord("test_del3", "chr1", 3000, true, "chr1", 3999, false, StructuralVariantType.DEL,
                1000, Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL),
                Collections.emptyList(), Collections.emptyMap());
        Assert.assertFalse(defragmenter.clusterTogether(deletion, deletion4), "Should barely not cluster");
    }

    @DataProvider(name = "testOverlapBuilders")
    public Object[][] testOverlapBuilders() {
        return new Object[][]{
            {
                new GenotypeBuilder()
                        .attribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT, 2)
                        .alleles(Lists.newArrayList(Allele.NO_CALL, Allele.NO_CALL)),
                    new GenotypeBuilder()
                            .attribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT, 1)
                            .alleles(Lists.newArrayList(Allele.NO_CALL, Allele.NO_CALL)),
                    new GenotypeBuilder()
                            .attribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT, 0)
                            .alleles(Lists.newArrayList(Allele.NO_CALL, Allele.NO_CALL)),
                    "CN only"
            },
            {
                new GenotypeBuilder()
                        .alleles(Lists.newArrayList(Allele.REF_N, Allele.REF_N)),
                    new GenotypeBuilder()
                            .alleles(Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL)),
                    new GenotypeBuilder()
                            .alleles(Lists.newArrayList(Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DEL)),
                    "GT only"
            }
        };
    }

    @Test(dataProvider= "testOverlapBuilders")
    private void testSampleOverlap(final GenotypeBuilder refGenotypeMaker,
                                   final GenotypeBuilder altGenotypeMaker1,
                                   final GenotypeBuilder altGenotypeMaker2,
                                   final String name) {
        // alt - ref - ref
        final SVCallRecord deletionARR1 = new SVCallRecord("test_del_arr1", "chr1", 1000, true, "chr1", 1999, false, StructuralVariantType.DEL,
                1000, Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL),
                Lists.newArrayList(altGenotypeMaker1.name("sample1").make(), refGenotypeMaker.name("sample2").make(), refGenotypeMaker.name("sample3").make()),
                Collections.emptyMap());
        // identical to above
        final SVCallRecord deletionARR1Copy = new SVCallRecord("test_del_arr1_copy", "chr1", 1000, true, "chr1", 1999, false, StructuralVariantType.DEL,
                1000, Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL),
                Lists.newArrayList(altGenotypeMaker1.name("sample1").make(), refGenotypeMaker.name("sample2").make(), refGenotypeMaker.name("sample3").make()),
                Collections.emptyMap());
        // second alt genotype - single-sample case will not cluster with ARR1
        final SVCallRecord deletionARR2 = new SVCallRecord("test_del_arr2", "chr1", 1000, true, "chr1", 1999, false, StructuralVariantType.DEL,
                1000, Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL),
                Lists.newArrayList(altGenotypeMaker2.name("sample1").make(), refGenotypeMaker.name("sample2").make(), refGenotypeMaker.name("sample3").make()),
                Collections.emptyMap());

        // ref - alt - ref
        final SVCallRecord deletionRAR = new SVCallRecord("test_del_rar", "chr1", 1000, true, "chr1", 1999, false, StructuralVariantType.DEL,
                1000, Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL),
                Lists.newArrayList(refGenotypeMaker.name("sample1").make(), altGenotypeMaker1.name("sample2").make(), refGenotypeMaker.name("sample3").make()),
                Collections.emptyMap());
        // alt - alt - ref
        final SVCallRecord deletionAAR = new SVCallRecord("test_del_aar", "chr1", 1000, true, "chr1", 1999, false, StructuralVariantType.DEL,
                1000, Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL),
                Lists.newArrayList(altGenotypeMaker1.name("sample1").make(), altGenotypeMaker1.name("sample2").make(), refGenotypeMaker.name("sample3").make()),
                Collections.emptyMap());
        // identical to above but second alt genotype - will cluster because there are multiple samples
        final SVCallRecord deletionAAR2 = new SVCallRecord("test_del_aar", "chr1", 1000, true, "chr1", 1999, false, StructuralVariantType.DEL,
                1000, Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL),
                Lists.newArrayList(altGenotypeMaker1.name("sample1").make(), altGenotypeMaker1.name("sample2").make(), refGenotypeMaker.name("sample3").make()),
                Collections.emptyMap());

        // alt - alt - alt
        final SVCallRecord deletionAAA = new SVCallRecord("test_del_aaa", "chr1", 1000, true, "chr1", 1999, false, StructuralVariantType.DEL,
                1000, Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL),
                Lists.newArrayList(altGenotypeMaker1.name("sample1").make(), altGenotypeMaker1.name("sample2").make(), altGenotypeMaker1.name("sample3").make()),
                Collections.emptyMap());
        Assert.assertFalse(defragmenter.clusterTogether(deletionARR1, deletionARR2));
        Assert.assertFalse(defragmenter.clusterTogether(deletionARR1, deletionRAR));
        Assert.assertFalse(defragmenter.clusterTogether(deletionARR1, deletionAAR));
        Assert.assertFalse(defragmenter.clusterTogether(deletionARR1, deletionAAA));
        Assert.assertFalse(defragmenter.clusterTogether(deletionRAR, deletionAAR));
        Assert.assertFalse(defragmenter.clusterTogether(deletionRAR, deletionAAA));

        Assert.assertTrue(defragmenter.clusterTogether(deletionARR1, deletionARR1Copy));
        Assert.assertTrue(defragmenter.clusterTogether(deletionAAR, deletionAAA));
        Assert.assertTrue(defragmenter.clusterTogether(deletionAAR, deletionAAR2));
    }

    @DataProvider(name = "maxPositionIntervals")
    public Object[][] recordPairs() {
        return new Object[][]{
                {100, 200},
                {50000, 500000},
                {1, 1},
                {1, 2},
        };
    }

    @Test(dataProvider= "maxPositionIntervals")
    public void testGetMaxClusterableStartingPosition(final int start, final int end) {
        final SVCallRecord call1 = new SVCallRecord("call1", "chr1", start, true, "chr1", end, false, StructuralVariantType.DEL,
                end - start + 1, Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL),
                Collections.emptyList(), Collections.emptyMap());
        final int maxClusterableStart = defragmenter.getMaxClusterableStartingPosition(call1);

        final int call2Start = maxClusterableStart;
        final int call2End = dictionary.getSequence(call1.getContigA()).getSequenceLength();
        final SVCallRecord call2 = new SVCallRecord("call2", "chr1", call2Start, true, "chr1", call2End, false, StructuralVariantType.DEL,
                call2End - call2Start + 1, Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL),
                Collections.emptyList(), Collections.emptyMap());
        Assert.assertTrue(defragmenter.clusterTogether(call1, call2));

        final int call3Start = maxClusterableStart + 1;
        final int call3End = dictionary.getSequence(call1.getContigA()).getSequenceLength();
        final SVCallRecord call3 = new SVCallRecord("call3", "chr1", call3Start, true, "chr1", call3End, false, StructuralVariantType.DEL,
                call3End - call3Start + 1, Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL),
                Collections.emptyList(), Collections.emptyMap());
        Assert.assertFalse(defragmenter.clusterTogether(call1, call3));
    }

}