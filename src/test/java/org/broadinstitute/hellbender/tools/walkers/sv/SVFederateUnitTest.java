package org.broadinstitute.hellbender.tools.walkers.sv;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SVTestUtils;
import org.broadinstitute.hellbender.tools.sv.cluster.CanonicalSVCollapser;
import org.broadinstitute.hellbender.tools.sv.cluster.CanonicalSVLinkage;
import org.broadinstitute.hellbender.tools.sv.cluster.SVFederationCollapser;
import org.testng.annotations.BeforeMethod;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class SVFederateUnitTest {

    private static final String SOURCE_A = "sourceA";
    private static final String SOURCE_B = "sourceB";
    private SVFederate federate;

    @BeforeMethod
    public void setUp() {
        federate = new SVFederate();
        federate.dictionary = SVTestUtils.hg38Dict;
        federate.sexes = new ArrayList<>(List.of("XY"));
        federate.afGroupingsA = new ArrayList<>(List.of("EUR"));
        federate.afGroupingsB = new ArrayList<>(List.of("AFR"));
        federate.afGroupingsAll = new ArrayList<>(List.of("EUR", "AFR"));
        federate.sourceToPrefixMap = new HashMap<>();
        federate.sourceToPrefixMap.put(SOURCE_A, "A");
        federate.sourceToPrefixMap.put(SOURCE_B, "B");
        federate.sourceToAFGroupingsMap = new HashMap<>();
        federate.sourceToAFGroupingsMap.put(SOURCE_A, List.of("EUR"));
        federate.sourceToAFGroupingsMap.put(SOURCE_B, List.of("AFR"));
        federate.collapser = new SVFederationCollapser(
                SVTestUtils.hg38Reference,
                CanonicalSVCollapser.AltAlleleSummaryStrategy.MOST_SPECIFIC_SUBTYPE,
                CanonicalSVCollapser.BreakpointSummaryStrategy.REPRESENTATIVE,
                CanonicalSVCollapser.FlagFieldLogic.OR);
        federate.linkage = new CanonicalSVLinkage<>(SVTestUtils.hg38Dict, true);
    }

    @Test
    public void testComputeLogAlleleFrequencyDifference() {
        Assert.assertEquals(SVFederate.computeLogAlleleFrequencyDifference(0.1, 0.01), Double.valueOf(1.0));
        Assert.assertTrue(Double.isNaN(SVFederate.computeLogAlleleFrequencyDifference(Double.NaN, 0.1)));
    }

    @Test
    public void testGetInfoKeyWithGroups() {
        Assert.assertEquals(federate.getInfoKeyWithGroups(VCFConstants.ALLELE_COUNT_KEY, "", ""), VCFConstants.ALLELE_COUNT_KEY);
        Assert.assertEquals(federate.getInfoKeyWithGroups(VCFConstants.ALLELE_COUNT_KEY, "EUR", ""), VCFConstants.ALLELE_COUNT_KEY + "_EUR");
        Assert.assertEquals(federate.getInfoKeyWithGroups(VCFConstants.ALLELE_COUNT_KEY, "EUR", "XY"), VCFConstants.ALLELE_COUNT_KEY + "_EUR_XY");
    }

    @Test
    public void testAnnotateUnmergedVariantForBiallelicAndCnv() {
        final Map<String, Object> biallelicAttributes = new HashMap<>();
        biallelicAttributes.put(GATKSVVCFConstants.SVTYPE, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL);
        biallelicAttributes.put(GATKSVVCFConstants.ALGORITHMS_ATTRIBUTE, SVTestUtils.PESR_ONLY_ALGORITHM_LIST);
        biallelicAttributes.put(VCFConstants.ALLELE_FREQUENCY_KEY, 0.2);
        biallelicAttributes.put(VCFConstants.ALLELE_COUNT_KEY, 4);
        biallelicAttributes.put(VCFConstants.ALLELE_NUMBER_KEY, 20);
        biallelicAttributes.put("N_HET", 2);
        biallelicAttributes.put("N_HOMREF", 10);
        biallelicAttributes.put("N_HOMALT", 0);
        biallelicAttributes.put(VCFConstants.ALLELE_FREQUENCY_KEY + "_XY", 0.2);
        biallelicAttributes.put(VCFConstants.ALLELE_COUNT_KEY + "_XY", 1);
        biallelicAttributes.put(VCFConstants.ALLELE_NUMBER_KEY + "_XY", 5);
        biallelicAttributes.put("N_HET_XY", 0);
        biallelicAttributes.put("N_HOMREF_XY", 4);
        biallelicAttributes.put("N_HOMALT_XY", 0);
        biallelicAttributes.put("N_HEMIALT_XY", 1);
        biallelicAttributes.put(VCFConstants.ALLELE_FREQUENCY_KEY + "_EUR", 0.2);
        biallelicAttributes.put(VCFConstants.ALLELE_COUNT_KEY + "_EUR", 4);
        biallelicAttributes.put(VCFConstants.ALLELE_NUMBER_KEY + "_EUR", 20);
        biallelicAttributes.put("N_HET_EUR", 2);
        biallelicAttributes.put("N_HOMREF_EUR", 10);
        biallelicAttributes.put("N_HOMALT_EUR", 0);
        biallelicAttributes.put(VCFConstants.ALLELE_FREQUENCY_KEY + "_EUR_XY", 0.2);
        biallelicAttributes.put(VCFConstants.ALLELE_COUNT_KEY + "_EUR_XY", 1);
        biallelicAttributes.put(VCFConstants.ALLELE_NUMBER_KEY + "_EUR_XY", 5);
        biallelicAttributes.put("N_HET_EUR_XY", 0);
        biallelicAttributes.put("N_HOMREF_EUR_XY", 4);
        biallelicAttributes.put("N_HOMALT_EUR_XY", 0);
        biallelicAttributes.put("N_HEMIALT_EUR_XY", 1);
        final SVCallRecord biallelic = federate.annotateUnmergedVariant(new VariantContextBuilder()
                .source(SOURCE_A)
                .id("biallelicA")
                .chr("chr1")
                .start(1000)
                .stop(1999)
                .alleles(List.of(Allele.REF_N, Allele.SV_SIMPLE_DEL))
                .attributes(biallelicAttributes)
                .make());
        Assert.assertEquals(biallelic.getAttributes().get("A_VID"), "biallelicA");
        Assert.assertEquals(biallelic.getAttributes().get("A_SVTYPE"), GATKSVVCFConstants.StructuralVariantAnnotationType.DEL);
        Assert.assertEquals(biallelic.getAttributes().get("A_AF"), 0.2);
        Assert.assertEquals(biallelic.getAttributes().get("A_AF_XY"), 0.2);
        Assert.assertEquals(biallelic.getAttributes().get("A_AF_EUR"), 0.2);
        Assert.assertEquals(biallelic.getAttributes().get("A_AF_EUR_XY"), 0.2);
        Assert.assertFalse(biallelic.getAttributes().containsKey("A_CN_NUMBER"));

        final Map<String, Object> cnvAttributes = new HashMap<>();
        cnvAttributes.put(GATKSVVCFConstants.SVTYPE, GATKSVVCFConstants.StructuralVariantAnnotationType.CNV);
        cnvAttributes.put(GATKSVVCFConstants.ALGORITHMS_ATTRIBUTE, SVTestUtils.DEPTH_ONLY_ALGORITHM_LIST);
        cnvAttributes.put("CN_NUMBER", 6);
        cnvAttributes.put("CN_COUNT", List.of(2, 3, 1));
        cnvAttributes.put("CN_NONREF_COUNT", 4);
        cnvAttributes.put("CN_NONREF_FREQ", 4.0 / 6.0);
        cnvAttributes.put("CN_FREQ", new double[]{2.0 / 6.0, 3.0 / 6.0, 1.0 / 6.0});
        cnvAttributes.put("CN_NUMBER_AFR", 6);
        cnvAttributes.put("CN_COUNT_AFR", List.of(2, 3, 1));
        cnvAttributes.put("CN_NONREF_COUNT_AFR", 4);
        cnvAttributes.put("CN_NONREF_FREQ_AFR", 4.0 / 6.0);
        cnvAttributes.put("CN_FREQ_AFR", new double[]{2.0 / 6.0, 3.0 / 6.0, 1.0 / 6.0});
        final SVCallRecord cnv = federate.annotateUnmergedVariant(new VariantContextBuilder()
            .source(SOURCE_B)
            .id("cnvB")
            .chr("chr1")
            .start(1000)
            .stop(1999)
            .alleles(List.of(Allele.REF_N, Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DUP))
            .attributes(cnvAttributes)
            .make());
        Assert.assertEquals(cnv.getAttributes().get("B_VID"), "cnvB");
        Assert.assertEquals(cnv.getAttributes().get("B_SVTYPE"), GATKSVVCFConstants.StructuralVariantAnnotationType.CNV);
        Assert.assertEquals(cnv.getAttributes().get("B_CN_NUMBER"), 6);
        Assert.assertEquals(cnv.getAttributes().get("B_CN_NUMBER_AFR"), 6);
        Assert.assertFalse(cnv.getAttributes().containsKey("B_AF"));
    }

    @Test
    public void testMergeBiallelicVariants() {
    final Map<String, Object> leftAttributes = new HashMap<>();
    leftAttributes.put(GATKSVVCFConstants.SVTYPE, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL);
    leftAttributes.put(GATKSVVCFConstants.ALGORITHMS_ATTRIBUTE, SVTestUtils.PESR_ONLY_ALGORITHM_LIST);
    leftAttributes.put(VCFConstants.ALLELE_FREQUENCY_KEY, 0.2);
    leftAttributes.put(VCFConstants.ALLELE_COUNT_KEY, 4);
    leftAttributes.put(VCFConstants.ALLELE_NUMBER_KEY, 20);
    leftAttributes.put("N_HET", 2);
    leftAttributes.put("N_HOMREF", 10);
    leftAttributes.put("N_HOMALT", 0);
    leftAttributes.put(VCFConstants.ALLELE_FREQUENCY_KEY + "_XY", 0.2);
    leftAttributes.put(VCFConstants.ALLELE_COUNT_KEY + "_XY", 1);
    leftAttributes.put(VCFConstants.ALLELE_NUMBER_KEY + "_XY", 5);
    leftAttributes.put("N_HET_XY", 0);
    leftAttributes.put("N_HOMREF_XY", 4);
    leftAttributes.put("N_HOMALT_XY", 0);
    leftAttributes.put("N_HEMIALT_XY", 1);
    leftAttributes.put(VCFConstants.ALLELE_FREQUENCY_KEY + "_EUR", 0.2);
    leftAttributes.put(VCFConstants.ALLELE_COUNT_KEY + "_EUR", 4);
    leftAttributes.put(VCFConstants.ALLELE_NUMBER_KEY + "_EUR", 20);
    leftAttributes.put("N_HET_EUR", 2);
    leftAttributes.put("N_HOMREF_EUR", 10);
    leftAttributes.put("N_HOMALT_EUR", 0);
    leftAttributes.put(VCFConstants.ALLELE_FREQUENCY_KEY + "_EUR_XY", 0.2);
    leftAttributes.put(VCFConstants.ALLELE_COUNT_KEY + "_EUR_XY", 1);
    leftAttributes.put(VCFConstants.ALLELE_NUMBER_KEY + "_EUR_XY", 5);
    leftAttributes.put("N_HET_EUR_XY", 0);
    leftAttributes.put("N_HOMREF_EUR_XY", 4);
    leftAttributes.put("N_HOMALT_EUR_XY", 0);
    leftAttributes.put("N_HEMIALT_EUR_XY", 1);
    final VariantContext left = new VariantContextBuilder()
        .source(SOURCE_A)
        .id("left")
        .chr("chr1")
        .start(1000)
        .stop(1999)
        .alleles(List.of(Allele.REF_N, Allele.SV_SIMPLE_DEL))
        .attributes(leftAttributes)
        .make();

    final Map<String, Object> rightAttributes = new HashMap<>();
    rightAttributes.put(GATKSVVCFConstants.SVTYPE, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL);
    rightAttributes.put(GATKSVVCFConstants.ALGORITHMS_ATTRIBUTE, SVTestUtils.PESR_ONLY_ALGORITHM_LIST);
    rightAttributes.put(VCFConstants.ALLELE_FREQUENCY_KEY, 0.1);
    rightAttributes.put(VCFConstants.ALLELE_COUNT_KEY, 2);
    rightAttributes.put(VCFConstants.ALLELE_NUMBER_KEY, 20);
    rightAttributes.put("N_HET", 1);
    rightAttributes.put("N_HOMREF", 10);
    rightAttributes.put("N_HOMALT", 0);
    rightAttributes.put(VCFConstants.ALLELE_FREQUENCY_KEY + "_XY", 0.1);
    rightAttributes.put(VCFConstants.ALLELE_COUNT_KEY + "_XY", 1);
    rightAttributes.put(VCFConstants.ALLELE_NUMBER_KEY + "_XY", 5);
    rightAttributes.put("N_HET_XY", 0);
    rightAttributes.put("N_HOMREF_XY", 4);
    rightAttributes.put("N_HOMALT_XY", 0);
    rightAttributes.put("N_HEMIALT_XY", 1);
    rightAttributes.put(VCFConstants.ALLELE_FREQUENCY_KEY + "_AFR", 0.1);
    rightAttributes.put(VCFConstants.ALLELE_COUNT_KEY + "_AFR", 2);
    rightAttributes.put(VCFConstants.ALLELE_NUMBER_KEY + "_AFR", 20);
    rightAttributes.put("N_HET_AFR", 1);
    rightAttributes.put("N_HOMREF_AFR", 10);
    rightAttributes.put("N_HOMALT_AFR", 0);
    rightAttributes.put(VCFConstants.ALLELE_FREQUENCY_KEY + "_AFR_XY", 0.1);
    rightAttributes.put(VCFConstants.ALLELE_COUNT_KEY + "_AFR_XY", 1);
    rightAttributes.put(VCFConstants.ALLELE_NUMBER_KEY + "_AFR_XY", 5);
    rightAttributes.put("N_HET_AFR_XY", 0);
    rightAttributes.put("N_HOMREF_AFR_XY", 4);
    rightAttributes.put("N_HOMALT_AFR_XY", 0);
    rightAttributes.put("N_HEMIALT_AFR_XY", 1);
    final VariantContext right = new VariantContextBuilder()
        .source(SOURCE_B)
        .id("right")
        .chr("chr1")
        .start(1000)
        .stop(1999)
        .alleles(List.of(Allele.REF_N, Allele.SV_SIMPLE_DEL))
        .attributes(rightAttributes)
        .make();

        final SVCallRecord merged = federate.merge(left, right);
        final Map<String, Object> attributes = merged.getAttributes();

        Assert.assertEquals(attributes.get("A_VID"), "left");
        Assert.assertEquals(attributes.get("B_VID"), "right");
        Assert.assertEquals(attributes.get(VCFConstants.ALLELE_COUNT_KEY), 6);
        Assert.assertEquals(attributes.get(VCFConstants.ALLELE_NUMBER_KEY), 40);
        Assert.assertEquals(attributes.get(VCFConstants.ALLELE_FREQUENCY_KEY), 0.15);
        Assert.assertEquals(attributes.get("N_HEMIALT_XY"), 2);
        Assert.assertEquals(attributes.get("AC_EUR"), 4);
        Assert.assertEquals(attributes.get("AC_AFR"), 2);
        Assert.assertEquals(attributes.get("LOG_AF_DIFFERENCE"), Double.valueOf(Math.abs(Math.log10(0.2) - Math.log10(0.1))));
        Assert.assertNotNull(attributes.get(GATKSVVCFConstants.RECIPROCAL_OVERLAP_INFO));
        Assert.assertNotNull(attributes.get(GATKSVVCFConstants.SIZE_SIMILARITY_INFO));
    }

    @Test
    public void testMergeMultiallelicVariantSetsBiallelicFieldsMissing() {
    final Map<String, Object> leftAttributes = new HashMap<>();
    leftAttributes.put(GATKSVVCFConstants.SVTYPE, GATKSVVCFConstants.StructuralVariantAnnotationType.CNV);
    leftAttributes.put(GATKSVVCFConstants.ALGORITHMS_ATTRIBUTE, SVTestUtils.DEPTH_ONLY_ALGORITHM_LIST);
    leftAttributes.put("CN_NUMBER", 6);
    leftAttributes.put("CN_COUNT", List.of(2, 3, 1));
    leftAttributes.put("CN_NONREF_COUNT", 4);
    leftAttributes.put("CN_NONREF_FREQ", 4.0 / 6.0);
    leftAttributes.put("CN_FREQ", new double[]{2.0 / 6.0, 3.0 / 6.0, 1.0 / 6.0});
    leftAttributes.put("CN_NUMBER_EUR", 6);
    leftAttributes.put("CN_COUNT_EUR", List.of(2, 3, 1));
    leftAttributes.put("CN_NONREF_COUNT_EUR", 4);
    leftAttributes.put("CN_NONREF_FREQ_EUR", 4.0 / 6.0);
    leftAttributes.put("CN_FREQ_EUR", new double[]{2.0 / 6.0, 3.0 / 6.0, 1.0 / 6.0});
    final VariantContext left = new VariantContextBuilder()
        .source(SOURCE_A)
        .id("cnvA")
        .chr("chr1")
        .start(1000)
        .stop(1999)
        .alleles(List.of(Allele.REF_N, Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DUP))
        .attributes(leftAttributes)
        .make();

    final Map<String, Object> rightAttributes = new HashMap<>();
    rightAttributes.put(GATKSVVCFConstants.SVTYPE, GATKSVVCFConstants.StructuralVariantAnnotationType.CNV);
    rightAttributes.put(GATKSVVCFConstants.ALGORITHMS_ATTRIBUTE, SVTestUtils.DEPTH_ONLY_ALGORITHM_LIST);
    rightAttributes.put("CN_NUMBER", 8);
    rightAttributes.put("CN_COUNT", List.of(1, 2, 5));
    rightAttributes.put("CN_NONREF_COUNT", 7);
    rightAttributes.put("CN_NONREF_FREQ", 7.0 / 8.0);
    rightAttributes.put("CN_FREQ", new double[]{1.0 / 8.0, 2.0 / 8.0, 5.0 / 8.0});
    rightAttributes.put("CN_NUMBER_AFR", 8);
    rightAttributes.put("CN_COUNT_AFR", List.of(1, 2, 5));
    rightAttributes.put("CN_NONREF_COUNT_AFR", 7);
    rightAttributes.put("CN_NONREF_FREQ_AFR", 7.0 / 8.0);
    rightAttributes.put("CN_FREQ_AFR", new double[]{1.0 / 8.0, 2.0 / 8.0, 5.0 / 8.0});
    final VariantContext right = new VariantContextBuilder()
        .source(SOURCE_B)
        .id("cnvB")
        .chr("chr1")
        .start(1000)
        .stop(1999)
        .alleles(List.of(Allele.REF_N, Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DUP))
        .attributes(rightAttributes)
        .make();

        final SVCallRecord merged = federate.merge(left, right);
        final Map<String, Object> attributes = merged.getAttributes();

        Assert.assertEquals(attributes.get(VCFConstants.ALLELE_COUNT_KEY), VCFConstants.MISSING_VALUE_v4);
        Assert.assertEquals(attributes.get("LOG_AF_DIFFERENCE"), VCFConstants.MISSING_VALUE_v4);
        Assert.assertEquals(attributes.get("A_CN_NUMBER"), 6);
        Assert.assertEquals(attributes.get("B_CN_NUMBER"), 8);
        Assert.assertEquals(attributes.get("CN_NUMBER"), 14);

        final int[] counts = (int[]) attributes.get("CN_COUNT");
        Assert.assertEquals(counts, new int[]{3, 5, 6});
        final double[] freqs = (double[]) attributes.get("CN_FREQ");
        Assert.assertEquals(freqs.length, 3);
        Assert.assertEquals(freqs[0], 3.0 / 14.0, 1e-10);
    }

    @Test
    public void testBuildVariantContextRemovesClusterMemberIdsAndClearsGenotypes() {
        federate.variantPrefix = "fed";

        final Map<String, Object> attributes = new HashMap<>();
        attributes.put(GATKSVVCFConstants.CLUSTER_MEMBER_IDS_KEY, List.of("member1", "member2"));
        final SVCallRecord record = new SVCallRecord("original", "chr1", 1000, true, "chr1", 1999, false,
                GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, null, Collections.emptyList(), 1000,
                Collections.emptyList(), SVTestUtils.DEPTH_ONLY_ALGORITHM_LIST,
                List.of(Allele.REF_N, Allele.SV_SIMPLE_DEL), SVTestUtils.threeGenotypes,
                attributes, Collections.emptySet(), null, SVTestUtils.hg38Dict);

        final VariantContext variant = federate.buildVariantContext(record);

        Assert.assertEquals(variant.getID(), "fed00000000");
        Assert.assertFalse(variant.hasAttribute(GATKSVVCFConstants.CLUSTER_MEMBER_IDS_KEY));
        Assert.assertTrue(variant.getGenotypes().isEmpty());
    }
}