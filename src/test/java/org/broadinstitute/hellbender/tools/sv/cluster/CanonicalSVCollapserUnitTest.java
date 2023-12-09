package org.broadinstitute.hellbender.tools.sv.cluster;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SVCallRecordUtils;
import org.broadinstitute.hellbender.tools.sv.SVTestUtils;
import org.broadinstitute.hellbender.utils.variant.GATKSVVariantContextUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class CanonicalSVCollapserUnitTest {

    private static final CanonicalSVCollapser collapser = new CanonicalSVCollapser(
            SVTestUtils.hg38Reference,
            CanonicalSVCollapser.AltAlleleSummaryStrategy.COMMON_SUBTYPE,
            CanonicalSVCollapser.BreakpointSummaryStrategy.MEDIAN_START_MEDIAN_END);
    private static final CanonicalSVCollapser collapserMinMax = new CanonicalSVCollapser(
            SVTestUtils.hg38Reference,
            CanonicalSVCollapser.AltAlleleSummaryStrategy.COMMON_SUBTYPE,
            CanonicalSVCollapser.BreakpointSummaryStrategy.MIN_START_MAX_END);
    private static final CanonicalSVCollapser collapserMaxMin = new CanonicalSVCollapser(
            SVTestUtils.hg38Reference,
            CanonicalSVCollapser.AltAlleleSummaryStrategy.COMMON_SUBTYPE,
            CanonicalSVCollapser.BreakpointSummaryStrategy.MAX_START_MIN_END);
    private static final CanonicalSVCollapser collapserMean = new CanonicalSVCollapser(
            SVTestUtils.hg38Reference,
            CanonicalSVCollapser.AltAlleleSummaryStrategy.COMMON_SUBTYPE,
            CanonicalSVCollapser.BreakpointSummaryStrategy.MEAN_START_MEAN_END);
    private static final CanonicalSVCollapser collapserRepresentative = new CanonicalSVCollapser(
            SVTestUtils.hg38Reference,
            CanonicalSVCollapser.AltAlleleSummaryStrategy.COMMON_SUBTYPE,
            CanonicalSVCollapser.BreakpointSummaryStrategy.REPRESENTATIVE);
    private static final CanonicalSVCollapser collapserSpecificAltAllele = new CanonicalSVCollapser(
            SVTestUtils.hg38Reference,
            CanonicalSVCollapser.AltAlleleSummaryStrategy.MOST_SPECIFIC_SUBTYPE,
            CanonicalSVCollapser.BreakpointSummaryStrategy.MEDIAN_START_MEDIAN_END);

    private static final Allele MEI_INSERTION_ALLELE = Allele.create("<INS:MEI>");
    private static final Allele SVA_INSERTION_ALLELE = Allele.create("<INS:MEI:SVA>");
    private static final Allele LINE_INSERTION_ALLELE = Allele.create("<INS:MEI:LINE>");

    @DataProvider(name = "clusterData")
    public Object[][] getClusterData() {
        return new Object[][]{
                // One deletion
                {
                        new SVClusterEngine.OutputCluster(
                                Lists.newArrayList(
                                        SVTestUtils.makeRecord("record1", "chr1", 1000, true,
                                                "chr1", 2000, false, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL,
                                                null, Collections.emptyList(), Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL),
                                                Lists.newArrayList(
                                                        new GenotypeBuilder("sample1", Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL)).attribute(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, 2),
                                                        new GenotypeBuilder("sample2", Lists.newArrayList(Allele.REF_N, Allele.REF_N)).attribute(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, 2)
                                                )
                                        )
                                )
                        ),
                        SVTestUtils.makeRecord("record1", "chr1", 1000, true,
                                "chr1", 2000, false, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL,
                                null, Collections.emptyList(), Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL),
                                Lists.newArrayList(
                                        new GenotypeBuilder("sample1", Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL)).attribute(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, 2),
                                        new GenotypeBuilder("sample2", Lists.newArrayList(Allele.REF_N, Allele.REF_N)).attribute(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, 2)
                                )
                        )
                },
                // Two deletions
                {
                    new SVClusterEngine.OutputCluster(
                        Lists.newArrayList(
                                SVTestUtils.makeRecord("record1", "chr1", 1000, true,
                                        "chr1", 2000, false, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL,
                                        null, Collections.emptyList(), Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL),
                                        Lists.newArrayList(
                                                new GenotypeBuilder("sample1", Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL)).attribute(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, 2),
                                                new GenotypeBuilder("sample2", Lists.newArrayList(Allele.REF_N, Allele.REF_N)).attribute(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, 2)
                                        )
                                ),
                                SVTestUtils.makeRecord("record2", "chr1", 1000, true,
                                        "chr1", 2000, false, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL,
                                        null, Collections.emptyList(), Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL),
                                        Lists.newArrayList(
                                                new GenotypeBuilder("sample1", Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL)).attribute(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, 2),
                                                new GenotypeBuilder("sample2", Lists.newArrayList(Allele.REF_N, Allele.REF_N)).attribute(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, 2)
                                        )
                                )
                        )
                    ),
                    SVTestUtils.makeRecord("record1", "chr1", 1000, true,
                            "chr1", 2000, false, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL,
                            null, Collections.emptyList(), Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL),
                            Lists.newArrayList(
                                    new GenotypeBuilder("sample1", Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL)).attribute(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, 2),
                                    new GenotypeBuilder("sample2", Lists.newArrayList(Allele.REF_N, Allele.REF_N)).attribute(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, 2)
                            )
                    )
                },
                // Breakends
                {
                        new SVClusterEngine.OutputCluster(
                                Lists.newArrayList(
                                        SVTestUtils.makeRecord("record1", "chr1", 1000, true,
                                                "chr2", 2000, false, GATKSVVCFConstants.StructuralVariantAnnotationType.BND,
                                                null, Collections.emptyList(), Lists.newArrayList(Allele.REF_N, GATKSVVariantContextUtils.BND_ALLELE),
                                                Lists.newArrayList(
                                                        new GenotypeBuilder("sample1", Lists.newArrayList(Allele.REF_N, GATKSVVariantContextUtils.BND_ALLELE)).attribute(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, 2),
                                                        new GenotypeBuilder("sample2", Lists.newArrayList(Allele.REF_N, Allele.REF_N)).attribute(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, 2)
                                                )
                                        ),
                                        SVTestUtils.makeRecord("record2", "chr1", 1000, true,
                                                "chr2", 2000, false, GATKSVVCFConstants.StructuralVariantAnnotationType.BND,
                                                null, Collections.emptyList(), Lists.newArrayList(Allele.REF_N, GATKSVVariantContextUtils.BND_ALLELE),
                                                Lists.newArrayList(
                                                        new GenotypeBuilder("sample1", Lists.newArrayList(Allele.REF_N, GATKSVVariantContextUtils.BND_ALLELE)).attribute(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, 2),
                                                        new GenotypeBuilder("sample2", Lists.newArrayList(Allele.REF_N, Allele.REF_N)).attribute(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, 2)
                                                )
                                        )
                                )
                        ),
                        SVTestUtils.makeRecord("record1", "chr1", 1000, true,
                                "chr2", 2000, false, GATKSVVCFConstants.StructuralVariantAnnotationType.BND,
                                null, Collections.emptyList(), Lists.newArrayList(Allele.REF_N, Allele.create("N[chr2:2000[")),
                                Lists.newArrayList(
                                        new GenotypeBuilder("sample1", Lists.newArrayList(Allele.REF_N, Allele.create("N[chr2:2000["))).attribute(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, 2),
                                        new GenotypeBuilder("sample2", Lists.newArrayList(Allele.REF_N, Allele.REF_N)).attribute(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, 2)
                                )
                        )
                }
        };
    }

    @Test(dataProvider = "clusterData")
    public void testCollapse(SVClusterEngine.OutputCluster cluster, SVCallRecord expected) {
        SVCallRecord result = collapser.collapse(cluster);
        SVTestUtils.assertEqualsExceptExcludedAttributes(result, expected, Collections.singletonList(GATKSVVCFConstants.CLUSTER_MEMBER_IDS_KEY));
    }

    @DataProvider(name = "collapseRefAllelesTestData")
    public Object[][] collapseRefAllelesTestData() {
        return new Object[][]{
                {1, Allele.REF_N},
                {100000, Allele.REF_C},
                {100001, Allele.REF_A},
                {100002, Allele.REF_C},
                {100003, Allele.REF_T},
                {248956422, Allele.REF_N},
        };
    }

    @Test(dataProvider= "collapseRefAllelesTestData")
    public void testCollapseRefAlleles(final int pos, final Allele result) {
        Assert.assertEquals(collapser.collapseRefAlleles("chr1", pos), result);
    }

    @DataProvider(name = "harmonizeAltAllelesTestData")
    public Object[][] harmonizeAltAllelesTestData() {
        return new Object[][]{
                // Basic cases
                {Collections.singletonList(Allele.SV_SIMPLE_DEL),
                        Collections.singletonList(new GenotypeBuilder().alleles(Lists.newArrayList(Allele.REF_A, Allele.SV_SIMPLE_DEL)).make()),
                        null},
                {Collections.singletonList(Allele.SV_SIMPLE_DEL),
                        Collections.singletonList(new GenotypeBuilder().alleles(Lists.newArrayList(Allele.REF_A, Allele.REF_A)).make()),
                        null},
                {Collections.emptyList(),
                        Collections.singletonList(new GenotypeBuilder().alleles(Lists.newArrayList(Allele.REF_A, Allele.REF_A)).make()),
                        null},
                // Multiallelic, without subtypes
                {Lists.newArrayList(Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DUP),
                        Collections.singletonList(new GenotypeBuilder().alleles(Lists.newArrayList(Allele.REF_A, Allele.SV_SIMPLE_DEL)).make()),
                        null},
                {Lists.newArrayList(Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DUP),
                        Collections.singletonList(new GenotypeBuilder().alleles(Lists.newArrayList(Allele.REF_A, Allele.REF_A)).make()),
                        null},
                {Lists.newArrayList(Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DUP),
                        Lists.newArrayList(
                                new GenotypeBuilder().alleles(Lists.newArrayList(Allele.REF_A, Allele.SV_SIMPLE_DEL)).make(),
                                new GenotypeBuilder().alleles(Lists.newArrayList(Allele.REF_A, Allele.SV_SIMPLE_DUP)).make()
                        ),
                        null},
                // Biallelic, with subtypes
                {Collections.singletonList(Allele.SV_SIMPLE_DUP),
                        Collections.singletonList(new GenotypeBuilder().alleles(Lists.newArrayList(Allele.REF_A, Allele.create("<DUP:TANDEM>"))).make()),
                        Collections.singletonList(new GenotypeBuilder().alleles(Lists.newArrayList(Allele.REF_A, Allele.SV_SIMPLE_DUP)).make())},
                {Collections.singletonList(Allele.create("<DUP:TANDEM>")),
                        Collections.singletonList(new GenotypeBuilder().alleles(Lists.newArrayList(Allele.REF_A, Allele.SV_SIMPLE_DUP)).make()),
                        Collections.singletonList(new GenotypeBuilder().alleles(Lists.newArrayList(Allele.REF_A, Allele.create("<DUP:TANDEM>"))).make())},
                {Collections.singletonList(Allele.create("<DUP:TANDEM>")),
                        Lists.newArrayList(
                                new GenotypeBuilder().alleles(Lists.newArrayList(Allele.REF_A, Allele.SV_SIMPLE_DUP)).make(),
                                new GenotypeBuilder().alleles(Lists.newArrayList(Allele.REF_A, Allele.create("<DUP:TANDEM>"))).make()
                        ),
                        Lists.newArrayList(
                                new GenotypeBuilder().alleles(Lists.newArrayList(Allele.REF_A, Allele.create("<DUP:TANDEM>"))).make(),
                                new GenotypeBuilder().alleles(Lists.newArrayList(Allele.REF_A, Allele.create("<DUP:TANDEM>"))).make()
                        )},
                {Lists.newArrayList(Allele.create("<INS:MEI>")),
                        Lists.newArrayList(
                                new GenotypeBuilder().alleles(Lists.newArrayList(Allele.REF_A, Allele.SV_SIMPLE_INS)).make(),
                                new GenotypeBuilder().alleles(Lists.newArrayList(Allele.REF_A, Allele.create("<INS:MEI>"))).make()
                        ),
                        Lists.newArrayList(
                                new GenotypeBuilder().alleles(Lists.newArrayList(Allele.REF_A, Allele.create("<INS:MEI>"))).make(),
                                new GenotypeBuilder().alleles(Lists.newArrayList(Allele.REF_A, Allele.create("<INS:MEI>"))).make()
                        )},
        };
    }

    @Test(dataProvider= "harmonizeAltAllelesTestData")
    public void testHarmonizeAltAlleles(final List<Allele> altAlleles, final List<Genotype> genotypes,
                                        final List<Genotype> expectedOrNull) {
        final List<Genotype> result = collapser.harmonizeAltAlleles(altAlleles, genotypes);
        // If null input, just expect the original genotypes
        final List<Genotype> expected = expectedOrNull == null ? genotypes : expectedOrNull;
        Assert.assertEquals(result.size(), expected.size());
        for (int i = 0; i < expected.size(); i++) {
            VariantContextTestUtils.assertGenotypesAreEqual(result.get(i), expected.get(i));
        }
    }

    @DataProvider(name = "collapseAltAllelesTestData")
    public Object[][] collapseAltAllelesTestData() {
        return new Object[][]{
                {
                    Collections.singletonList(Collections.emptyList()),
                        GATKSVVCFConstants.StructuralVariantAnnotationType.DEL,
                        new Integer[]{0},
                        new Integer[]{0},
                        Collections.emptyList(),
                        Collections.emptyList()
                },
                {
                    Collections.singletonList(Collections.singletonList(Allele.REF_N)),
                        GATKSVVCFConstants.StructuralVariantAnnotationType.DEL,
                        new Integer[]{1},
                        new Integer[]{1},
                        Collections.emptyList(),
                        Collections.emptyList()
                },
                {
                    Collections.singletonList(Collections.singletonList(Allele.REF_A)),
                        GATKSVVCFConstants.StructuralVariantAnnotationType.DEL,
                        new Integer[]{1},
                        new Integer[]{1},
                        Collections.emptyList(),
                        Collections.emptyList()
                },
                {
                    Collections.singletonList(Collections.singletonList(Allele.NO_CALL)),
                        GATKSVVCFConstants.StructuralVariantAnnotationType.DEL,
                        new Integer[]{1},
                        new Integer[]{1},
                        Collections.emptyList(),
                        Collections.emptyList()
                },
                {
                    Collections.singletonList(Collections.singletonList(Allele.ALT_A)),
                        GATKSVVCFConstants.StructuralVariantAnnotationType.INS,
                        new Integer[]{1},
                        new Integer[]{null},
                        Collections.singletonList(Allele.ALT_A),
                        Collections.singletonList(Allele.ALT_A)
                },
                {
                    Collections.singletonList(Collections.singletonList(Allele.SV_SIMPLE_DEL)),
                        GATKSVVCFConstants.StructuralVariantAnnotationType.DEL,
                        new Integer[]{1},
                        new Integer[]{0},
                        Collections.singletonList(Allele.SV_SIMPLE_DEL),
                        Collections.singletonList(Allele.SV_SIMPLE_DEL)
                },
                {
                    Lists.newArrayList(Collections.singletonList(Allele.REF_N), Collections.singletonList(Allele.SV_SIMPLE_DEL)),
                        GATKSVVCFConstants.StructuralVariantAnnotationType.DEL,
                        new Integer[]{1, 1},
                        new Integer[]{1, 0},
                        Collections.singletonList(Allele.SV_SIMPLE_DEL),
                        Collections.singletonList(Allele.SV_SIMPLE_DEL)
                },
                {
                    Lists.newArrayList(Lists.newArrayList(Allele.REF_N, Allele.REF_N), Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL)),
                        GATKSVVCFConstants.StructuralVariantAnnotationType.DEL,
                        new Integer[]{2, 2},
                        new Integer[]{2, 1},
                        Collections.singletonList(Allele.SV_SIMPLE_DEL),
                        Collections.singletonList(Allele.SV_SIMPLE_DEL)
                },
                {
                    Lists.newArrayList(Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DUP), Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL)),
                        GATKSVVCFConstants.StructuralVariantAnnotationType.CNV,
                        new Integer[]{2, 2},
                        new Integer[]{3, 1},
                        Lists.newArrayList(Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DUP),
                        Lists.newArrayList(Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DUP)},
                {
                    Lists.newArrayList(Collections.singletonList(MEI_INSERTION_ALLELE), Collections.singletonList(MEI_INSERTION_ALLELE)),
                        GATKSVVCFConstants.StructuralVariantAnnotationType.INS,
                        new Integer[]{1, 1},
                        new Integer[]{null, null},
                        Collections.singletonList(MEI_INSERTION_ALLELE),
                        Collections.singletonList(MEI_INSERTION_ALLELE)},
                {
                    Lists.newArrayList(Collections.singletonList(Allele.SV_SIMPLE_INS), Collections.singletonList(MEI_INSERTION_ALLELE)),
                        GATKSVVCFConstants.StructuralVariantAnnotationType.INS,
                        new Integer[]{1, 1},
                        new Integer[]{null, null},
                        Collections.singletonList(Allele.SV_SIMPLE_INS),
                        Collections.singletonList(MEI_INSERTION_ALLELE)},
                {
                    Lists.newArrayList(Collections.singletonList(MEI_INSERTION_ALLELE), Collections.singletonList(SVA_INSERTION_ALLELE)),
                        GATKSVVCFConstants.StructuralVariantAnnotationType.INS,
                        new Integer[]{1, 1},
                        new Integer[]{null, null},
                        Collections.singletonList(MEI_INSERTION_ALLELE),
                        Collections.singletonList(SVA_INSERTION_ALLELE)},
                {
                    Lists.newArrayList(Collections.singletonList(LINE_INSERTION_ALLELE), Collections.singletonList(SVA_INSERTION_ALLELE)),
                        GATKSVVCFConstants.StructuralVariantAnnotationType.INS,
                        new Integer[]{1, 1},
                        new Integer[]{null, null},
                        Collections.singletonList(MEI_INSERTION_ALLELE),
                        Collections.singletonList(MEI_INSERTION_ALLELE)},
        };
    }

    @Test(dataProvider= "collapseAltAllelesTestData")
    public void collapseAltAllelesTest(final List<List<Allele>> recordGenotypeAlleles,
                                       final GATKSVVCFConstants.StructuralVariantAnnotationType svtype,
                                       final Integer[] expectedCopyNumber,
                                       final Integer[] copyNumber,
                                       final List<Allele> resultCommon,
                                       final List<Allele> resultSpecific) {
        final List<Allele> variantAlleles = recordGenotypeAlleles.stream().flatMap(List::stream).distinct().collect(Collectors.toList());
        final List<SVCallRecord> records = IntStream.range(0, recordGenotypeAlleles.size())
                .mapToObj(i -> SVTestUtils.newCallRecordWithAlleles(recordGenotypeAlleles.get(i), variantAlleles, svtype, expectedCopyNumber[i], copyNumber[i]))
                .collect(Collectors.toList());

        final List<Allele> sortedTestCommon = SVCallRecordUtils.sortAlleles(collapser.collapseAltAlleles(records));
        final List<Allele> sortedExpectedCommon = SVCallRecordUtils.sortAlleles(resultCommon);
        Assert.assertEquals(sortedTestCommon, sortedExpectedCommon);

        final List<Allele> sortedTestSpecific = SVCallRecordUtils.sortAlleles(collapserSpecificAltAllele.collapseAltAlleles(records));
        final List<Allele> sortedExpectedSpecific = SVCallRecordUtils.sortAlleles(resultSpecific);
        Assert.assertEquals(sortedTestSpecific, sortedExpectedSpecific);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void collapseInvalidAltAllelesTest() {
        final List<SVCallRecord> records = Lists.newArrayList(
                SVTestUtils.newCallRecordWithAlleles(
                        Lists.newArrayList(Allele.REF_N, Allele.create("N[chr1:1[", false)),
                        Lists.newArrayList(Allele.REF_N, Allele.create("N[chr1:1[", false)),
                        GATKSVVCFConstants.StructuralVariantAnnotationType.INS,
                        2, 2
                ),
                SVTestUtils.newCallRecordWithAlleles(
                        Lists.newArrayList(Allele.REF_N, Allele.create("N[chr1:2[", false)),
                        Lists.newArrayList(Allele.REF_N, Allele.create("N[chr1:2[", false)),
                        GATKSVVCFConstants.StructuralVariantAnnotationType.INS,
                        2, 2
                )
        );
        collapser.collapseAltAlleles(records);
    }

    @DataProvider(name = "bndAlleleData")
    public Object[][] bndAlleleData() {
        return new Object[][] {
                { true, true, "contigB", 10, "A", "A]contigB:10]" },
                { true, false, "contigB", 10, "T", "T[contigB:10[" },
                { false, true, "contigB", 20, "C", "]contigB:20]C" },
                { false, false, "contigB", 20, "G", "[contigB:20[G" },
        };
    }

    @Test(dataProvider = "bndAlleleData")
    public void testConstructBndAllele(Boolean strandA, Boolean strandB, String contigB, int posB,
                                       String refAlleleString, String expectedAlleleString) {
        final Allele refAllele = Allele.create(refAlleleString, true);
        final Allele expected = Allele.create(expectedAlleleString, false);
        Allele result = CanonicalSVCollapser.constructBndAllele(strandA, strandB, contigB, posB, refAllele);
        Assert.assertEquals(result, expected);
    }

    private static final String TEST_KEY_1 = "TEST_KEY_1";

    private Map<String, Object> createGenotypeTestAttributes(final Integer expectedCopyNumber, final String testVal) {
        final String[] keys = new String[]{
                GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT,
                TEST_KEY_1
        };
        final Object[] vals = new Object[]{expectedCopyNumber, testVal};
        return SVTestUtils.keyValueArraysToMap(keys, vals);
    }

    private Map<String, Object> createGenotypeTestAttributes(final Integer expectedCopyNumber, final Integer copyNumber) {
        final String[] keys = new String[]{
                GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT,
                GATKSVVCFConstants.COPY_NUMBER_FORMAT
        };
        final Object[] vals = new Object[]{expectedCopyNumber, copyNumber};
        return SVTestUtils.keyValueArraysToMap(keys, vals);
    }

    private Map<String, Object> createGenotypeTestAttributesWithGQ(final Integer expectedCopyNumber, final Integer genotypeQuality) {
        final String[] keys = new String[]{
                GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT,
                VCFConstants.GENOTYPE_QUALITY_KEY
        };
        final Object[] vals = new Object[]{expectedCopyNumber, genotypeQuality};
        return SVTestUtils.keyValueArraysToMap(keys, vals);
    }

    private Map<String, Object> createGenotypeTestAttributesWithCNQ(final Integer expectedCopyNumber, final Integer copyNumber, final Integer copyNumberQuality) {
        final String[] keys = new String[]{
                GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT,
                GATKSVVCFConstants.COPY_NUMBER_FORMAT,
                GATKSVVCFConstants.COPY_NUMBER_QUALITY_FORMAT
        };
        final Object[] vals = new Object[]{expectedCopyNumber, copyNumber, copyNumberQuality};
        return SVTestUtils.keyValueArraysToMap(keys, vals);
    }

    private Map<String, Object> createGenotypeTestAttributes(final Integer expectedCopyNumber) {
        return Collections.singletonMap(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, expectedCopyNumber);
    }

    @DataProvider(name = "collapseSampleGenotypesTestData")
    public Object[][] collapseSampleGenotypesTestData() {
        return new Object[][]{
                // Empty case
                {
                        "sample",
                        Collections.singletonList(Collections.emptyList()),
                        Collections.singletonList(createGenotypeTestAttributes(0)),
                        Allele.REF_N,
                        GenotypeBuilder.create(
                                "sample",
                                Collections.emptyList(),
                                createGenotypeTestAttributes(0)
                        )
                },
                // Extra attribute
                {
                        "sample",
                        Collections.singletonList(Collections.emptyList()),
                        Collections.singletonList(createGenotypeTestAttributes(0, "test")),
                        Allele.REF_N,
                        GenotypeBuilder.create(
                                "sample",
                                Collections.emptyList(),
                                createGenotypeTestAttributes(0, "test")
                        )
                },
                // Haploid no-call
                {
                        "sample",
                        Lists.newArrayList(
                                Collections.singletonList(Allele.NO_CALL),
                                Collections.singletonList(Allele.NO_CALL)
                        ),
                        Lists.newArrayList(
                                createGenotypeTestAttributes(1),
                                createGenotypeTestAttributes(1)
                        ),
                        Allele.REF_N,
                        GenotypeBuilder.create(
                                "sample",
                                Collections.singletonList(Allele.NO_CALL),
                                createGenotypeTestAttributes(1)
                        )
                },
                // Simple ref haploid, different ref allele
                {
                        "sample",
                        Lists.newArrayList(
                                Collections.singletonList(Allele.REF_T),
                                Collections.singletonList(Allele.REF_T)
                        ),
                        Lists.newArrayList(
                                createGenotypeTestAttributes(1),
                                createGenotypeTestAttributes(1)
                        ),
                        Allele.REF_T,
                        GenotypeBuilder.create(
                                "sample",
                                Collections.singletonList(Allele.REF_T),
                                createGenotypeTestAttributes(1)
                        )
                },
                // Simple ref haploid, with no-call
                {
                        "sample",
                        Lists.newArrayList(
                                Collections.singletonList(Allele.REF_T),
                                Collections.singletonList(Allele.NO_CALL)
                        ),
                        Lists.newArrayList(
                                createGenotypeTestAttributes(1),
                                createGenotypeTestAttributes(1)
                        ),
                        Allele.REF_T,
                        GenotypeBuilder.create(
                                "sample",
                                Collections.singletonList(Allele.REF_T),
                                createGenotypeTestAttributes(1)
                        )
                },
                // Simple ref diploid
                {
                        "sample",
                        Lists.newArrayList(
                                Lists.newArrayList(Allele.REF_N, Allele.REF_N),
                                Lists.newArrayList(Allele.REF_N, Allele.REF_N)
                        ),
                        Lists.newArrayList(
                                createGenotypeTestAttributes(2),
                                createGenotypeTestAttributes(2)
                        ),
                        Allele.REF_N,
                        GenotypeBuilder.create(
                                "sample",
                                Lists.newArrayList(Allele.REF_N, Allele.REF_N),
                                createGenotypeTestAttributes(2)
                        )
                },
                // Simple ref diploid, with no-call
                {
                        "sample",
                        Lists.newArrayList(
                                Lists.newArrayList(Allele.REF_N, Allele.REF_N),
                                Lists.newArrayList(Allele.NO_CALL, Allele.NO_CALL)
                        ),
                        Lists.newArrayList(
                                createGenotypeTestAttributes(2),
                                createGenotypeTestAttributes(2)
                        ),
                        Allele.REF_N,
                        GenotypeBuilder.create(
                                "sample",
                                Lists.newArrayList(Allele.REF_N, Allele.REF_N),
                                createGenotypeTestAttributes(2)
                        )
                },
                // Simple INS cases
                {
                        "sample",
                        Lists.newArrayList(
                                Collections.singletonList(Allele.REF_N),
                                Collections.singletonList(Allele.SV_SIMPLE_INS)
                        ),
                        Lists.newArrayList(
                                createGenotypeTestAttributes(1),
                                createGenotypeTestAttributes(1)
                        ),
                        Allele.REF_N,
                        GenotypeBuilder.create(
                                "sample",
                                Lists.newArrayList(Allele.SV_SIMPLE_INS),
                                createGenotypeTestAttributes(1)
                        )
                },
                // Prefer non-ref over higher GQ
                {
                        "sample",
                        Lists.newArrayList(
                                Collections.singletonList(Allele.NO_CALL),
                                Collections.singletonList(Allele.SV_SIMPLE_INS)
                        ),
                        Lists.newArrayList(
                                createGenotypeTestAttributesWithGQ(1, 30),
                                createGenotypeTestAttributesWithGQ(1, 20)
                        ),
                        Allele.REF_N,
                        GenotypeBuilder.create(
                                "sample",
                                Lists.newArrayList(Allele.SV_SIMPLE_INS),
                                createGenotypeTestAttributesWithGQ(1, 20)
                        )
                },
                // Use higher GQ
                {
                        "sample",
                        Lists.newArrayList(
                                Collections.singletonList(Allele.SV_SIMPLE_INS),
                                Collections.singletonList(Allele.SV_SIMPLE_INS)
                        ),
                        Lists.newArrayList(
                                createGenotypeTestAttributesWithGQ(1, 20),
                                createGenotypeTestAttributesWithGQ(1, 30)
                        ),
                        Allele.REF_N,
                        GenotypeBuilder.create(
                                "sample",
                                Lists.newArrayList(Allele.SV_SIMPLE_INS),
                                createGenotypeTestAttributesWithGQ(1, 30)
                        )
                },
                // het preferred over hom ref
                {
                        "sample",
                        Lists.newArrayList(
                                Lists.newArrayList(Allele.REF_N, Allele.REF_N),
                                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_INS)
                        ),
                        Lists.newArrayList(
                                createGenotypeTestAttributes(2),
                                createGenotypeTestAttributes(2)
                        ),
                        Allele.REF_N,
                        GenotypeBuilder.create(
                                "sample",
                                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_INS),
                                createGenotypeTestAttributes(2)
                        )
                },
                // het preferred over hom ref even with lower gq
                {
                        "sample",
                        Lists.newArrayList(
                                Lists.newArrayList(Allele.REF_N, Allele.REF_N),
                                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_INS)
                        ),
                        Lists.newArrayList(
                                createGenotypeTestAttributesWithGQ(2, 30),
                                createGenotypeTestAttributesWithGQ(2, 20)
                        ),
                        Allele.REF_N,
                        GenotypeBuilder.create(
                                "sample",
                                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_INS),
                                createGenotypeTestAttributesWithGQ(2, 20)
                        )
                },
                // hom var preferred over het
                {
                        "sample",
                        Lists.newArrayList(
                                Lists.newArrayList(Allele.SV_SIMPLE_INS, Allele.SV_SIMPLE_INS),
                                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_INS)
                        ),
                        Lists.newArrayList(
                                createGenotypeTestAttributes(2),
                                createGenotypeTestAttributes(2)
                        ),
                        Allele.REF_N,
                        GenotypeBuilder.create(
                                "sample",
                                Lists.newArrayList(Allele.SV_SIMPLE_INS, Allele.SV_SIMPLE_INS),
                                createGenotypeTestAttributes(2)
                        )
                },
                // hom-var over het if GQ equal
                {
                        "sample",
                        Lists.newArrayList(
                                Lists.newArrayList(Allele.SV_SIMPLE_INS, Allele.SV_SIMPLE_INS),
                                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_INS)
                        ),
                        Lists.newArrayList(
                                createGenotypeTestAttributesWithGQ(2, 30),
                                createGenotypeTestAttributesWithGQ(2, 30)
                        ),
                        Allele.REF_N,
                        GenotypeBuilder.create(
                                "sample",
                                Lists.newArrayList(Allele.SV_SIMPLE_INS, Allele.SV_SIMPLE_INS),
                                createGenotypeTestAttributes(2)
                        )
                },
                // het over hom-var if GQ is higher
                {
                        "sample",
                        Lists.newArrayList(
                                Lists.newArrayList(Allele.SV_SIMPLE_INS, Allele.SV_SIMPLE_INS),
                                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_INS)
                        ),
                        Lists.newArrayList(
                                createGenotypeTestAttributesWithGQ(2, 30),
                                createGenotypeTestAttributesWithGQ(2, 40)
                        ),
                        Allele.REF_N,
                        GenotypeBuilder.create(
                                "sample",
                                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_INS),
                                createGenotypeTestAttributesWithGQ(2, 40)
                        )
                },
                // hom var preferred over hom-ref too
                {
                        "sample",
                        Lists.newArrayList(
                                Lists.newArrayList(Allele.SV_SIMPLE_INS, Allele.SV_SIMPLE_INS),
                                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_INS),
                                Lists.newArrayList(Allele.REF_N, Allele.REF_N)
                        ),
                        Lists.newArrayList(
                                createGenotypeTestAttributesWithGQ(2, 30),
                                createGenotypeTestAttributesWithGQ(2, 30),
                                createGenotypeTestAttributesWithGQ(2, 30)
                        ),
                        Allele.REF_N,
                        GenotypeBuilder.create(
                                "sample",
                                Lists.newArrayList(Allele.SV_SIMPLE_INS, Allele.SV_SIMPLE_INS),
                                createGenotypeTestAttributesWithGQ(2, 30)
                        )
                },
                // Take highest non-ref GQ
                {
                        "sample",
                        Lists.newArrayList(
                                Lists.newArrayList(Allele.SV_SIMPLE_INS, Allele.SV_SIMPLE_INS),
                                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_INS),
                                Lists.newArrayList(Allele.REF_N, Allele.REF_N),
                                Lists.newArrayList(Allele.NO_CALL, Allele.NO_CALL)
                        ),
                        Lists.newArrayList(
                                createGenotypeTestAttributesWithGQ(2, 30),
                                createGenotypeTestAttributesWithGQ(2, 40),
                                createGenotypeTestAttributesWithGQ(2, 50),
                                createGenotypeTestAttributesWithGQ(2, 50)
                        ),
                        Allele.REF_N,
                        GenotypeBuilder.create(
                                "sample",
                                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_INS),
                                createGenotypeTestAttributesWithGQ(2, 40)
                        )
                },
                // triploid - 1 alt
                {
                        "sample",
                        Lists.newArrayList(
                                Lists.newArrayList(Allele.REF_N, Allele.REF_N, Allele.SV_SIMPLE_INS),
                                Lists.newArrayList(Allele.REF_N, Allele.REF_N, Allele.REF_N)
                        ),
                        Lists.newArrayList(
                                createGenotypeTestAttributes(3),
                                createGenotypeTestAttributes(3)
                        ),
                        Allele.REF_N,
                        GenotypeBuilder.create(
                                "sample",
                                Lists.newArrayList(Allele.REF_N, Allele.REF_N, Allele.SV_SIMPLE_INS),
                                createGenotypeTestAttributes(3)
                        )
                },
                {
                        "sample",
                        Lists.newArrayList(
                                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_INS, Allele.SV_SIMPLE_INS),
                                Lists.newArrayList(Allele.REF_N, Allele.REF_N, Allele.SV_SIMPLE_INS),
                                Lists.newArrayList(Allele.REF_N, Allele.REF_N, Allele.REF_N)
                        ),
                        Lists.newArrayList(
                                createGenotypeTestAttributes(3),
                                createGenotypeTestAttributes(3),
                                createGenotypeTestAttributes(3)
                        ),
                        Allele.REF_N,
                        GenotypeBuilder.create(
                                "sample",
                                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_INS, Allele.SV_SIMPLE_INS),
                                createGenotypeTestAttributes(3)
                        )
                },
                // triploid - 2 alt
                {
                        "sample",
                        Lists.newArrayList(
                                Lists.newArrayList(Allele.REF_N, Allele.REF_N, Allele.SV_SIMPLE_INS),
                                Lists.newArrayList(Allele.REF_N, Allele.REF_N, Allele.REF_N),
                                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_INS, Allele.SV_SIMPLE_INS),
                                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_INS, Allele.SV_SIMPLE_INS)
                        ),
                        Lists.newArrayList(
                                createGenotypeTestAttributes(3),
                                createGenotypeTestAttributes(3),
                                createGenotypeTestAttributes(3),
                                createGenotypeTestAttributes(3)
                        ),
                        Allele.REF_N,
                        GenotypeBuilder.create(
                                "sample",
                                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_INS, Allele.SV_SIMPLE_INS),
                                createGenotypeTestAttributes(3)
                        )
                },
                // triploid - 3 alt
                {
                        "sample",
                        Lists.newArrayList(
                                Lists.newArrayList(Allele.SV_SIMPLE_INS, Allele.SV_SIMPLE_INS, Allele.SV_SIMPLE_INS),
                                Lists.newArrayList(Allele.SV_SIMPLE_INS, Allele.SV_SIMPLE_INS, Allele.SV_SIMPLE_INS),
                                Lists.newArrayList(Allele.SV_SIMPLE_INS, Allele.SV_SIMPLE_INS, Allele.SV_SIMPLE_INS),
                                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_INS, Allele.SV_SIMPLE_INS),
                                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_INS, Allele.SV_SIMPLE_INS),
                                Lists.newArrayList(Allele.REF_N, Allele.REF_N, Allele.SV_SIMPLE_INS),
                                Lists.newArrayList(Allele.REF_N, Allele.REF_N, Allele.REF_N)
                        ),
                        Lists.newArrayList(
                                createGenotypeTestAttributes(3),
                                createGenotypeTestAttributes(3),
                                createGenotypeTestAttributes(3),
                                createGenotypeTestAttributes(3),
                                createGenotypeTestAttributes(3),
                                createGenotypeTestAttributes(3),
                                createGenotypeTestAttributes(3)
                        ),
                        Allele.REF_N,
                        GenotypeBuilder.create(
                                "sample",
                                Lists.newArrayList(Allele.SV_SIMPLE_INS, Allele.SV_SIMPLE_INS, Allele.SV_SIMPLE_INS),
                                createGenotypeTestAttributes(3)
                        )
                },
                // Simple DEL
                {
                        "sample",
                        Lists.newArrayList(
                                Collections.singletonList(Allele.REF_N),
                                Collections.singletonList(Allele.SV_SIMPLE_DEL)
                        ),
                        Lists.newArrayList(
                                createGenotypeTestAttributes(1, 1),
                                createGenotypeTestAttributes(1, 0)
                        ),
                        Allele.REF_N,
                        GenotypeBuilder.create(
                                "sample",
                                Lists.newArrayList(Allele.SV_SIMPLE_DEL),
                                createGenotypeTestAttributes(1, 0)
                        )
                },
                // Simple DEL, prefer called genotype (despite CN indicating a possible event)
                {
                        "sample",
                        Lists.newArrayList(
                                Collections.singletonList(Allele.REF_N),
                                Collections.singletonList(Allele.NO_CALL)
                        ),
                        Lists.newArrayList(
                                createGenotypeTestAttributes(1, 1),
                                createGenotypeTestAttributes(1, 0)
                        ),
                        Allele.REF_N,
                        GenotypeBuilder.create(
                                "sample",
                                Lists.newArrayList(Allele.REF_N),
                                createGenotypeTestAttributes(1, 1)
                        )
                },
                // Simple DEL diploid het
                {
                        "sample",
                        Lists.newArrayList(
                                Lists.newArrayList(Allele.REF_N, Allele.REF_N),
                                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL)
                        ),
                        Lists.newArrayList(
                                createGenotypeTestAttributes(2, 2),
                                createGenotypeTestAttributes(2, 1)
                        ),
                        Allele.REF_N,
                        GenotypeBuilder.create(
                                "sample",
                                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL),
                                createGenotypeTestAttributes(2, 1)
                        )
                },
                // Simple DEL diploid hom var
                {
                        "sample",
                        Lists.newArrayList(
                                Lists.newArrayList(Allele.REF_N, Allele.REF_N),
                                Lists.newArrayList(Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DEL)
                        ),
                        Lists.newArrayList(
                                createGenotypeTestAttributes(2, 2),
                                createGenotypeTestAttributes(2, 0)
                        ),
                        Allele.REF_N,
                        GenotypeBuilder.create(
                                "sample",
                                Lists.newArrayList(Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DEL),
                                createGenotypeTestAttributes(2, 0)
                        )
                },

                // Simple DUP, haploid
                {
                        "sample",
                        Lists.newArrayList(
                                Collections.singletonList(Allele.REF_N),
                                Collections.singletonList(Allele.SV_SIMPLE_DUP)
                        ),
                        Lists.newArrayList(
                                createGenotypeTestAttributes(1, 1),
                                createGenotypeTestAttributes(1, 2)
                        ),
                        Allele.REF_N,
                        GenotypeBuilder.create(
                                "sample",
                                Lists.newArrayList(Allele.SV_SIMPLE_DUP),
                                createGenotypeTestAttributes(1, 2)
                        )
                },
                // Simple DUP, haploid, take higher CNQ
                {
                        "sample",
                        Lists.newArrayList(
                                Collections.singletonList(Allele.SV_SIMPLE_DUP),
                                Collections.singletonList(Allele.SV_SIMPLE_DUP)
                        ),
                        Lists.newArrayList(
                                createGenotypeTestAttributesWithCNQ(1, 2, 20),
                                createGenotypeTestAttributesWithCNQ(1, 2, 30)
                        ),
                        Allele.REF_N,
                        GenotypeBuilder.create(
                                "sample",
                                Lists.newArrayList(Allele.SV_SIMPLE_DUP),
                                createGenotypeTestAttributesWithCNQ(1, 2, 30)
                        )
                },
                // Simple DUP diploid het
                {
                        "sample",
                        Lists.newArrayList(
                                Lists.newArrayList(Allele.REF_N, Allele.REF_N),
                                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DUP)
                        ),
                        Lists.newArrayList(
                                createGenotypeTestAttributes(2, 2),
                                createGenotypeTestAttributes(2, 3)
                        ),
                        Allele.REF_N,
                        GenotypeBuilder.create(
                                "sample",
                                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DUP),
                                createGenotypeTestAttributes(2, 3)
                        )
                },
                // Simple DUP diploid hom var
                {
                        "sample",
                        Lists.newArrayList(
                                Lists.newArrayList(Allele.REF_N, Allele.REF_N),
                                Lists.newArrayList(Allele.SV_SIMPLE_DUP, Allele.SV_SIMPLE_DUP)
                        ),
                        Lists.newArrayList(
                                createGenotypeTestAttributes(2, 2),
                                createGenotypeTestAttributes(2, 4)
                        ),
                        Allele.REF_N,
                        GenotypeBuilder.create(
                                "sample",
                                Lists.newArrayList(Allele.SV_SIMPLE_DUP, Allele.SV_SIMPLE_DUP),
                                createGenotypeTestAttributes(2, 4)
                        )
                },
                // Simple DUP diploid hom ref
                {
                        "sample",
                        Lists.newArrayList(
                                Lists.newArrayList(Allele.REF_N, Allele.REF_N),
                                Lists.newArrayList(Allele.REF_N, Allele.REF_N)
                        ),
                        Lists.newArrayList(
                                createGenotypeTestAttributes(2, 2),
                                createGenotypeTestAttributes(2, 2)
                        ),
                        Allele.REF_N,
                        GenotypeBuilder.create(
                                "sample",
                                Lists.newArrayList(Allele.REF_N, Allele.REF_N),
                                createGenotypeTestAttributes(2, 2)
                        )
                },
                // Simple DUP triploid with 1 alt - unambiguous alleles
                {
                        "sample",
                        Lists.newArrayList(
                                Lists.newArrayList(Allele.REF_N, Allele.REF_N, Allele.REF_N),
                                Lists.newArrayList(Allele.REF_N, Allele.REF_N, Allele.SV_SIMPLE_DUP)
                        ),
                        Lists.newArrayList(
                                createGenotypeTestAttributes(3, 3),
                                createGenotypeTestAttributes(3, 4)
                        ),
                        Allele.REF_N,
                        GenotypeBuilder.create(
                                "sample",
                                Lists.newArrayList(Allele.REF_N, Allele.REF_N, Allele.SV_SIMPLE_DUP),
                                createGenotypeTestAttributes(3, 4)
                        )
                },
                // Simple DUP triploid with 2 alts
                {
                        "sample",
                        Lists.newArrayList(
                                Lists.newArrayList(Allele.REF_N, Allele.REF_N, Allele.REF_N),
                                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DUP, Allele.SV_SIMPLE_DUP)
                        ),
                        Lists.newArrayList(
                                createGenotypeTestAttributes(3, 3),
                                createGenotypeTestAttributes(3, 5)
                        ),
                        Allele.REF_N,
                        GenotypeBuilder.create(
                                "sample",
                                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DUP, Allele.SV_SIMPLE_DUP),
                                createGenotypeTestAttributes(3, 5)
                        )
                },
                // Simple DUP triploid where het should be prioritized over no-call
                {
                        "sample",
                        Lists.newArrayList(
                                Lists.newArrayList(Allele.REF_N, Allele.REF_N, Allele.REF_N),
                                Lists.newArrayList(Allele.REF_N, Allele.REF_N, Allele.SV_SIMPLE_DUP),
                                Lists.newArrayList(Allele.NO_CALL, Allele.NO_CALL, Allele.NO_CALL)
                        ),
                        Lists.newArrayList(
                                createGenotypeTestAttributes(3, 3),
                                createGenotypeTestAttributes(3, 4),
                                createGenotypeTestAttributes(3, 5)
                        ),
                        Allele.REF_N,
                        GenotypeBuilder.create(
                                "sample",
                                Lists.newArrayList(Allele.REF_N, Allele.REF_N, Allele.SV_SIMPLE_DUP),
                                createGenotypeTestAttributes(3, 4)
                        )
                },

                // Multi-allelic CNV, haploid hom ref
                {
                        "sample",
                        Lists.newArrayList(
                                Collections.singletonList(Allele.NO_CALL),
                                Collections.singletonList(Allele.NO_CALL)
                        ),
                        Lists.newArrayList(
                                createGenotypeTestAttributes(1, 1),
                                createGenotypeTestAttributes(1, 1)
                        ),
                        Allele.REF_N,
                        GenotypeBuilder.create(
                                "sample",
                                Lists.newArrayList(Allele.NO_CALL),
                                createGenotypeTestAttributes(1, 1)
                        )
                },
                // Multi-allelic CNV, rare case with equal CNQ take CN!=ECN
                {
                        "sample",
                        Lists.newArrayList(
                                Collections.singletonList(Allele.NO_CALL),
                                Collections.singletonList(Allele.NO_CALL)
                        ),
                        Lists.newArrayList(
                                createGenotypeTestAttributesWithCNQ(1, 1, 30),
                                createGenotypeTestAttributesWithCNQ(1, 0, 30)
                        ),
                        Allele.REF_N,
                        GenotypeBuilder.create(
                                "sample",
                                Lists.newArrayList(Allele.NO_CALL),
                                createGenotypeTestAttributesWithCNQ(1, 0, 30)
                        )
                },
                // Multi-allelic CNV, haploid dup
                {
                        "sample",
                        Lists.newArrayList(
                                Collections.singletonList(Allele.NO_CALL),
                                Collections.singletonList(Allele.NO_CALL)
                        ),
                        Lists.newArrayList(
                                createGenotypeTestAttributes(1, 1),
                                createGenotypeTestAttributes(1, 2)
                        ),
                        Allele.REF_N,
                        GenotypeBuilder.create(
                                "sample",
                                Lists.newArrayList(Allele.NO_CALL),
                                createGenotypeTestAttributes(1, 2)
                        )
                },
                // Multi-allelic CNV, when CNQ equal use CN!=ECN
                {
                        "sample",
                        Lists.newArrayList(
                                Lists.newArrayList(Allele.NO_CALL, Allele.NO_CALL),
                                Lists.newArrayList(Allele.NO_CALL, Allele.NO_CALL)
                        ),
                        Lists.newArrayList(
                                createGenotypeTestAttributesWithCNQ(2, 2, 30),
                                createGenotypeTestAttributesWithCNQ(2, 1, 30)
                        ),
                        Allele.REF_N,
                        GenotypeBuilder.create(
                                "sample",
                                Lists.newArrayList(Allele.NO_CALL, Allele.NO_CALL),
                                createGenotypeTestAttributesWithCNQ(2, 1, 30)
                        )
                },
                // Multi-allelic CNV, when CNQ equal use CN!=ECN
                {
                        "sample",
                        Lists.newArrayList(
                                Lists.newArrayList(Allele.NO_CALL, Allele.NO_CALL),
                                Lists.newArrayList(Allele.NO_CALL, Allele.NO_CALL)
                        ),
                        Lists.newArrayList(
                                createGenotypeTestAttributesWithCNQ(2, 2, 30),
                                createGenotypeTestAttributesWithCNQ(2, 0, 30)
                        ),
                        Allele.REF_N,
                        GenotypeBuilder.create(
                                "sample",
                                Lists.newArrayList(Allele.NO_CALL, Allele.NO_CALL),
                                createGenotypeTestAttributesWithCNQ(2, 0, 30)
                        )
                },
                // Multi-allelic CNV, conflicting del and dup genotypes determined by CNQ
                {
                        "sample",
                        Lists.newArrayList(
                                Collections.singletonList(Allele.NO_CALL),
                                Collections.singletonList(Allele.NO_CALL)
                        ),
                        Lists.newArrayList(
                                createGenotypeTestAttributesWithCNQ(1, 2, 30),
                                createGenotypeTestAttributesWithCNQ(1, 0, 40)
                        ),
                        Allele.REF_N,
                        GenotypeBuilder.create(
                                "sample",
                                Lists.newArrayList(Allele.NO_CALL),
                                createGenotypeTestAttributesWithCNQ(1, 0, 40)
                        )
                },
                {
                        "sample",
                        Lists.newArrayList(
                                Collections.singletonList(Allele.NO_CALL),
                                Collections.singletonList(Allele.NO_CALL)
                        ),
                        Lists.newArrayList(
                                createGenotypeTestAttributesWithCNQ(1, 2, 50),
                                createGenotypeTestAttributesWithCNQ(1, 0, 40)
                        ),
                        Allele.REF_N,
                        GenotypeBuilder.create(
                                "sample",
                                Lists.newArrayList(Allele.NO_CALL),
                                createGenotypeTestAttributesWithCNQ(1, 2, 50)
                        )
                },
        };
    }

    @Test(dataProvider="collapseSampleGenotypesTestData")
    public void collapseSampleGenotypesTest(final String sampleId,
                                            final List<List<Allele>> alleles,
                                            final List<Map<String, Object>> attributes,
                                            final Allele refAllele,
                                            final Genotype expected) {
        final List<Genotype> genotypes = IntStream.range(0, alleles.size())
                .mapToObj(i -> GenotypeBuilder.create(sampleId, alleles.get(i), attributes.get(i)))
                .collect(Collectors.toList());
        final Genotype test = collapser.collapseSampleGenotypes(genotypes, refAllele);
        VariantContextTestUtils.assertGenotypesAreEqual(test, expected);
    }

    @DataProvider(name = "makeBiallelicListTestData")
    public Object[][] makeBiallelicListTestData() {
        return new Object[][]{
                // No alleles
                {Allele.SV_SIMPLE_DEL, Allele.REF_A, 0, 0, Collections.emptyList()},
                {Allele.SV_SIMPLE_DEL, Allele.REF_A, 0, 1, Collections.singletonList(Allele.REF_A)},
                {Allele.SV_SIMPLE_DEL, Allele.REF_A, 1, 1, Collections.singletonList(Allele.SV_SIMPLE_DEL)},
                {Allele.SV_SIMPLE_INS, Allele.REF_A, 1, 1, Collections.singletonList(Allele.SV_SIMPLE_INS)},
                {Allele.SV_SIMPLE_DEL, Allele.REF_A, 1, 2, Lists.newArrayList(Allele.REF_A, Allele.SV_SIMPLE_DEL)},
                {Allele.SV_SIMPLE_DEL, Allele.REF_A, 2, 2, Lists.newArrayList(Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DEL)},
                {Allele.SV_SIMPLE_DEL, Allele.REF_A, 2, 2, Lists.newArrayList(Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DEL)},
                {Allele.SV_SIMPLE_DEL, Allele.REF_A, 2, 3, Lists.newArrayList(Allele.REF_A, Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DEL)},
        };
    }

    @Test(dataProvider="makeBiallelicListTestData")
    public void makeBiallelicListTest(final Allele alt, final Allele ref, final int numAlt,
                                      final int ploidy, final List<Allele> result) {
        final List<Allele> test = CanonicalSVCollapser.makeBiallelicList(alt, ref, numAlt, ploidy);
        Assert.assertEquals(test.stream().sorted().collect(Collectors.toList()),
                result.stream().sorted().collect(Collectors.toList()));
    }

    @DataProvider(name = "collapseLengthTestData")
    public Object[][] collapseLengthTestData() {
        return new Object[][]{
                { SVTestUtils.newCallRecordInsertionWithLength(100), GATKSVVCFConstants.StructuralVariantAnnotationType.INS, 100 },
                { SVTestUtils.newCallRecordInsertionWithLength(null), GATKSVVCFConstants.StructuralVariantAnnotationType.INS, null },
                { SVTestUtils.newBndCallRecordWithStrands(true, false), GATKSVVCFConstants.StructuralVariantAnnotationType.BND, null },
                { SVTestUtils.newCtxCallRecord(), GATKSVVCFConstants.StructuralVariantAnnotationType.CTX, null },
                { SVTestUtils.newCpxCallRecordWithLength(50), GATKSVVCFConstants.StructuralVariantAnnotationType.CPX, 50 },
                { SVTestUtils.newCallRecordWithLengthAndType(100, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL), GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, 100 },
                { SVTestUtils.newCallRecordWithLengthAndType(200, GATKSVVCFConstants.StructuralVariantAnnotationType.DUP), GATKSVVCFConstants.StructuralVariantAnnotationType.DUP, 200 },
                { SVTestUtils.newCallRecordWithLengthAndType(300, GATKSVVCFConstants.StructuralVariantAnnotationType.DUP), GATKSVVCFConstants.StructuralVariantAnnotationType.CNV, 300 }
        };
    }

    @Test(dataProvider= "collapseLengthTestData")
    public void collapseLengthTest(final SVCallRecord record,
                                   GATKSVVCFConstants.StructuralVariantAnnotationType type,
                                   final Integer expectedLength) {
        Assert.assertEquals(collapser.collapseLength(record, type, record.getPositionA(), record.getPositionB()), expectedLength);
    }

    @DataProvider(name = "collapseIdsTestData")
    public Object[][] collapseIdsTestData() {
        return new Object[][]{
                {Collections.singletonList("var1"), "var1"},
                {Lists.newArrayList("var1", "var2"), "var1"},
                {Lists.newArrayList("var2", "var1"), "var1"},
        };
    }

    @DataProvider(name = "getMostPreciseCallsTestData")
    public Object[][] getMostPreciseCallsTestData() {
        return new Object[][]{
                {
                        new String[]{"depth1"},
                        Collections.singletonList(
                                Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM)
                        ),
                        Sets.newHashSet("depth1")
                },
                {
                        new String[]{"depth1", "depth2"},
                        Lists.newArrayList(
                                Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
                                Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM)
                        ),
                        Sets.newHashSet("depth1", "depth2")
                },
                {
                        new String[]{"pesr1"},
                        Collections.singletonList(
                                SVTestUtils.PESR_ONLY_ALGORITHM_LIST
                        ),
                        Sets.newHashSet("pesr1")
                },
                {
                        new String[]{"pesr1", "pesr2"},
                        Lists.newArrayList(
                                SVTestUtils.PESR_ONLY_ALGORITHM_LIST,
                                SVTestUtils.PESR_ONLY_ALGORITHM_LIST
                        ),
                        Sets.newHashSet("pesr1", "pesr2")
                },
                {
                        new String[]{"depth1", "pesr1"},
                        Lists.newArrayList(
                                Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
                                SVTestUtils.PESR_ONLY_ALGORITHM_LIST
                        ),
                        Sets.newHashSet("pesr1"),
                },
                {
                        new String[]{"mixed1", "depth1"},
                        Lists.newArrayList(
                                Lists.newArrayList(GATKSVVCFConstants.DEPTH_ALGORITHM, SVTestUtils.PESR_ALGORITHM),
                                Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM)
                        ),
                        Sets.newHashSet("mixed1"),
                },
                {
                        new String[]{"mixed1", "pesr1"},
                        Lists.newArrayList(
                                Lists.newArrayList(GATKSVVCFConstants.DEPTH_ALGORITHM, SVTestUtils.PESR_ALGORITHM),
                                SVTestUtils.PESR_ONLY_ALGORITHM_LIST
                        ),
                        Sets.newHashSet("mixed1", "pesr1"),
                },
        };
    }

    @Test(dataProvider= "getMostPreciseCallsTestData")
    public void getMostPreciseCallsTest(final String[] ids, final List<List<String>> algorithms, final Set<String> expectedIds) {
        final List<SVCallRecord> records = IntStream.range(0, ids.length).mapToObj(i -> SVTestUtils.newDeletionCallRecordWithIdAndAlgorithms(ids[i], algorithms.get(i))).collect(Collectors.toList());
        final List<String> resultIds = collapser.getRecordsWithMostPreciseBreakpoints(records).stream().map(SVCallRecord::getId).collect(Collectors.toList());
        Assert.assertEquals(resultIds.size(), expectedIds.size());
        Assert.assertTrue(expectedIds.containsAll(resultIds));
    }

    @DataProvider(name = "collapseIntervalTestData")
    public Object[][] collapseIntervalTestData() {
        return new Object[][]{
                // 1 variant
                {
                        new String[] {"chr1", "chr1"},
                        new int[]{1001}, // starts
                        new int[]{1100}, // ends
                        GATKSVVCFConstants.StructuralVariantAnnotationType.DEL,
                        new int[]{1001, 1100}, // median strategy
                        new int[]{1001, 1100}, // min-max strategy
                        new int[]{1001, 1100}, // max-min strategy
                        new int[]{1001, 1100}  // mean strategy
                },
                // 2 variants
                {
                        new String[] {"chr1", "chr1"},
                        new int[]{1001, 1011},
                        new int[]{1100, 1110},
                        GATKSVVCFConstants.StructuralVariantAnnotationType.DEL,
                        new int[]{1001, 1100},
                        new int[]{1001, 1110},
                        new int[]{1011, 1100},
                        new int[]{1006, 1105}
                },
                // 2 variants
                {
                        new String[] {"chr1", "chr1"},
                        new int[]{1001, 1011},
                        new int[]{1000, 1110},
                        GATKSVVCFConstants.StructuralVariantAnnotationType.DEL,
                        new int[]{1001, 1001},  // true median  second position is 1000 but not allowed since < 1001
                        new int[]{1001, 1110},
                        new int[]{1011, 1011},  // true min second position is 1000 but not allowed since < 1011
                        new int[]{1006, 1055}
                },
                // 3 variants
                {
                        new String[] {"chr1", "chr1"},
                        new int[]{1001, 1011, 1021},
                        new int[]{1100, 1110, 1121},
                        GATKSVVCFConstants.StructuralVariantAnnotationType.DUP,
                        new int[]{1011, 1110},
                        new int[]{1001, 1121},
                        new int[]{1021, 1100},
                        new int[]{1011, 1110}
                },
                // 3 variants, ends overlap starts
                {
                        new String[] {"chr1", "chr1"},
                        new int[]{1001, 1011, 1021},
                        new int[]{1031, 1011, 1021},
                        GATKSVVCFConstants.StructuralVariantAnnotationType.DUP,
                        new int[]{1011, 1021},
                        new int[]{1001, 1031},
                        new int[]{1021, 1021}, // min end before max start
                        new int[]{1011, 1021}
                },
                // BND, same contig
                {
                        new String[] {"chr1", "chr1"},
                        new int[]{1001, 1011, 1021},
                        new int[]{1100, 1110, 1121},
                        GATKSVVCFConstants.StructuralVariantAnnotationType.BND,
                        new int[]{1011, 1110},
                        new int[]{1001, 1121},
                        new int[]{1021, 1100},
                        new int[]{1011, 1110}
                },
                // BND, different contigs
                {
                        new String[] {"chr1", "chr2"},
                        new int[]{1001, 1011, 1021},
                        new int[]{1000, 1110, 1121},
                        GATKSVVCFConstants.StructuralVariantAnnotationType.BND,
                        new int[]{1011, 1110},
                        new int[]{1001, 1121},
                        new int[]{1021, 1000}, // 1000 < 1021, but ok because on diff contigs
                        new int[]{1011, 1077}
                },
                // INS, same start/end
                {
                        new String[] {"chr1", "chr1"},
                        new int[]{1001, 1011, 1021},
                        new int[]{1001, 1011, 1021},
                        GATKSVVCFConstants.StructuralVariantAnnotationType.INS,
                        new int[]{1011, 1011},
                        new int[]{1001, 1001},
                        new int[]{1021, 1021},
                        new int[]{1011, 1011}
                },
                // INS, different start/end
                {
                        new String[] {"chr1", "chr1"},
                        new int[]{1001, 1011, 1021},
                        new int[]{1011, 1021, 1031},
                        GATKSVVCFConstants.StructuralVariantAnnotationType.INS,
                        new int[]{1011, 1011},
                        new int[]{1001, 1001},
                        new int[]{1021, 1021},
                        new int[]{1011, 1011}
                }
        };
    }

    @Test(dataProvider= "collapseIntervalTestData")
    public void collapseIntervalTest(final String[] contigs, final int[] starts, final int[] ends, final GATKSVVCFConstants.StructuralVariantAnnotationType svtype,
                                     final int[] expectedMedian, final int[] expectedMinMax, final int[] expectedMaxMin,
                                     final int[] expectedMean) {
        final List<SVCallRecord> records =  IntStream.range(0, starts.length)
                .mapToObj(i -> SVTestUtils.newCallRecordWithContigsIntervalAndType(contigs[0], starts[i], contigs[1], ends[i], svtype)).collect(Collectors.toList());
        collapseIntervalTestHelper(collapser, svtype, contigs, records, expectedMedian);
        collapseIntervalTestHelper(collapserMinMax, svtype, contigs, records, expectedMinMax);
        collapseIntervalTestHelper(collapserMaxMin, svtype, contigs, records, expectedMaxMin);
        collapseIntervalTestHelper(collapserMean, svtype, contigs, records, expectedMean);
    }

    @Test
    public void collapseIntervalRepresentativeTest() {
        // Choose second record with more carriers
        final List<SVCallRecord> records =
                Lists.newArrayList(
                        SVTestUtils.makeRecord("record1", "chr1", 1000, true,
                                "chr1", 2000, false, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL,
                                null, Collections.emptyList(), Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL),
                                Lists.newArrayList(
                                        new GenotypeBuilder("sample1", Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL)).attribute(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, 2),
                                        new GenotypeBuilder("sample2", Lists.newArrayList(Allele.REF_N, Allele.REF_N)).attribute(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, 2)
                                )
                        ),
                        SVTestUtils.makeRecord("record2", "chr1", 1001, true,
                                "chr1", 2001, false, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL,
                                null, Collections.emptyList(), Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL),
                                Lists.newArrayList(
                                        new GenotypeBuilder("sample1", Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL)).attribute(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, 2),
                                        new GenotypeBuilder("sample2", Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL)).attribute(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, 2)
                                )
                        )
                );
        final Pair<Integer, Integer> result = collapserRepresentative.collapseInterval(records);
        Assert.assertEquals((int) result.getLeft(), 1001);
        Assert.assertEquals((int) result.getRight(), 2001);

        // record2 and record3 have the best carrier status, but choose second record which is closer to all others on average
        final List<SVCallRecord> records2 =
                Lists.newArrayList(
                        SVTestUtils.makeRecord("record1", "chr1", 1000, true,
                                "chr1", 2000, false, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL,
                                null, Collections.emptyList(), Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL),
                                Lists.newArrayList(
                                        new GenotypeBuilder("sample1", Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL)).attribute(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, 2),
                                        new GenotypeBuilder("sample2", Lists.newArrayList(Allele.REF_N, Allele.REF_N)).attribute(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, 2)
                                )
                        ),
                        SVTestUtils.makeRecord("record2", "chr1", 999, true,
                                "chr1", 2000, false, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL,
                                null, Collections.emptyList(), Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL),
                                Lists.newArrayList(
                                        new GenotypeBuilder("sample1", Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL)).attribute(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, 2),
                                        new GenotypeBuilder("sample2", Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL)).attribute(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, 2)
                                )
                        ),
                        SVTestUtils.makeRecord("record3", "chr1", 1005, true,
                                "chr1", 2000, false, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL,
                                null, Collections.emptyList(), Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL),
                                Lists.newArrayList(
                                        new GenotypeBuilder("sample1", Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL)).attribute(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, 2),
                                        new GenotypeBuilder("sample2", Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL)).attribute(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, 2)
                                )
                        )
                );
        final Pair<Integer, Integer> result2 = collapserRepresentative.collapseInterval(records2);
        Assert.assertEquals((int) result2.getLeft(), 999);
        Assert.assertEquals((int) result2.getRight(), 2000);
    }

    @DataProvider(name = "distanceDataProvider")
    public Object[][] distanceDataProvider() {
        return new Object[][]{
                {5, 10, new int[]{0}, new int[]{0}, 15},
                {5, 10, new int[]{0}, new int[]{10}, 5},
                {5, 10, new int[]{5}, new int[]{0}, 10},
                {5, 10, new int[]{5}, new int[]{10}, 0},
                {5, 10, new int[]{1, 3, 7}, new int[]{9, 12}, 11},
                {0, 0, new int[]{0, 0, 0}, new int[]{0, 0}, 0},
                {-5, -10, new int[]{-1, -3, -7}, new int[]{-9, -12}, 11},
        };
    }

    @Test(dataProvider = "distanceDataProvider")
    public void testGetDistance(int posA, int posB, int[] starts, int[] ends, long expectedDistance) {
        final long actualDistance = CanonicalSVCollapser.getDistance(posA, posB, starts, ends);
        Assert.assertEquals(actualDistance, expectedDistance);
    }

    private static void collapseIntervalTestHelper(final CanonicalSVCollapser collapser,
                                                   final GATKSVVCFConstants.StructuralVariantAnnotationType svtype,
                                                   final String[] contigs,
                                                   final List<SVCallRecord> records,
                                                   final int[] expected) {
        final Pair<Integer, Integer> result = collapser.collapseInterval(records);
        Assert.assertEquals((int) result.getKey(), expected[0]);
        Assert.assertEquals((int) result.getValue(), expected[1]);
        // Test that the coordinates validate when we try to create a new record
        SVTestUtils.newCallRecordWithContigsIntervalAndType(contigs[0], result.getKey(), contigs[1], result.getValue(), svtype);
    }

    @DataProvider(name = "collapseTypesTestData")
    public Object[][] collapseTypesTestData() {
        return new Object[][]{
                {Collections.singletonList(GATKSVVCFConstants.StructuralVariantAnnotationType.DEL), GATKSVVCFConstants.StructuralVariantAnnotationType.DEL},
                {Collections.singletonList(GATKSVVCFConstants.StructuralVariantAnnotationType.DUP), GATKSVVCFConstants.StructuralVariantAnnotationType.DUP},
                {Collections.singletonList(GATKSVVCFConstants.StructuralVariantAnnotationType.INS), GATKSVVCFConstants.StructuralVariantAnnotationType.INS},
                {Collections.singletonList(GATKSVVCFConstants.StructuralVariantAnnotationType.INV), GATKSVVCFConstants.StructuralVariantAnnotationType.INV},
                {Collections.singletonList(GATKSVVCFConstants.StructuralVariantAnnotationType.BND), GATKSVVCFConstants.StructuralVariantAnnotationType.BND},
                {Collections.singletonList(GATKSVVCFConstants.StructuralVariantAnnotationType.CNV), GATKSVVCFConstants.StructuralVariantAnnotationType.CNV},
                {Lists.newArrayList(GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, GATKSVVCFConstants.StructuralVariantAnnotationType.DUP), GATKSVVCFConstants.StructuralVariantAnnotationType.CNV},
                {Lists.newArrayList(GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL), GATKSVVCFConstants.StructuralVariantAnnotationType.DEL},
                {Lists.newArrayList(GATKSVVCFConstants.StructuralVariantAnnotationType.DUP, GATKSVVCFConstants.StructuralVariantAnnotationType.DUP), GATKSVVCFConstants.StructuralVariantAnnotationType.DUP},
                {Lists.newArrayList(GATKSVVCFConstants.StructuralVariantAnnotationType.INS, GATKSVVCFConstants.StructuralVariantAnnotationType.INS), GATKSVVCFConstants.StructuralVariantAnnotationType.INS},
                {Lists.newArrayList(GATKSVVCFConstants.StructuralVariantAnnotationType.INV, GATKSVVCFConstants.StructuralVariantAnnotationType.INV), GATKSVVCFConstants.StructuralVariantAnnotationType.INV},
                {Lists.newArrayList(GATKSVVCFConstants.StructuralVariantAnnotationType.BND, GATKSVVCFConstants.StructuralVariantAnnotationType.BND), GATKSVVCFConstants.StructuralVariantAnnotationType.BND},
                {Lists.newArrayList(GATKSVVCFConstants.StructuralVariantAnnotationType.CNV, GATKSVVCFConstants.StructuralVariantAnnotationType.CNV), GATKSVVCFConstants.StructuralVariantAnnotationType.CNV},
                {Lists.newArrayList(GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, GATKSVVCFConstants.StructuralVariantAnnotationType.DUP, GATKSVVCFConstants.StructuralVariantAnnotationType.CNV), GATKSVVCFConstants.StructuralVariantAnnotationType.CNV}
        };
    }

    @Test(dataProvider= "collapseTypesTestData")
    public void collapseTypesTest(final List<GATKSVVCFConstants.StructuralVariantAnnotationType> types, final GATKSVVCFConstants.StructuralVariantAnnotationType expectedResult) {
        final List<SVCallRecord> records = types.stream().map(t -> SVTestUtils.newPESRCallRecordWithIntervalAndType(1, 100, t)).collect(Collectors.toList());
        Assert.assertEquals(collapser.collapseTypes(records), expectedResult);
    }

    @DataProvider(name = "collapseAlgorithmsTestData")
    public Object[][] collapseAlgorithmsTestData() {
        return new Object[][]{
                // depth only
                {Collections.singletonList(Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM)), Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM)},
                {Lists.newArrayList(Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM), Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM)),
                        Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM)},

                // pesr only
                {Collections.singletonList(SVTestUtils.PESR_ONLY_ALGORITHM_LIST), SVTestUtils.PESR_ONLY_ALGORITHM_LIST},
                {Collections.singletonList(Lists.newArrayList("pesrA", "pesrB")), Lists.newArrayList("pesrA", "pesrB")},
                {Collections.singletonList(Lists.newArrayList("pesrB", "pesrA")), Lists.newArrayList("pesrA", "pesrB")},
                {Lists.newArrayList(Collections.singletonList("pesrA"), Collections.singletonList("pesrB")), Lists.newArrayList("pesrA", "pesrB")},
                {Lists.newArrayList(Collections.singletonList("pesrB"), Collections.singletonList("pesrA")), Lists.newArrayList("pesrA", "pesrB")},
                {Lists.newArrayList(Lists.newArrayList("pesrA", "pesrB"), Collections.singletonList("pesrB")), Lists.newArrayList("pesrA", "pesrB")},
                {Lists.newArrayList(Lists.newArrayList("pesrB", "pesrA"), Collections.singletonList("pesrA")), Lists.newArrayList("pesrA", "pesrB")},

                // mixed depth/pesr
                {Lists.newArrayList(SVTestUtils.PESR_ONLY_ALGORITHM_LIST, SVTestUtils.DEPTH_ONLY_ALGORITHM_LIST),
                        Lists.newArrayList(GATKSVVCFConstants.DEPTH_ALGORITHM, SVTestUtils.PESR_ALGORITHM)},
                {Lists.newArrayList(Lists.newArrayList(GATKSVVCFConstants.DEPTH_ALGORITHM, SVTestUtils.PESR_ALGORITHM),
                        Lists.newArrayList(GATKSVVCFConstants.DEPTH_ALGORITHM, SVTestUtils.PESR_ALGORITHM)),
                        Lists.newArrayList(GATKSVVCFConstants.DEPTH_ALGORITHM, SVTestUtils.PESR_ALGORITHM)}
        };
    }

    @Test(dataProvider= "collapseAlgorithmsTestData")
    public void collapseAlgorithmsTest(final List<List<String>> algorithmLists, final List<String> expectedResult) {
        final List<SVCallRecord> records = algorithmLists.stream().map(list -> SVTestUtils.newDeletionCallRecordWithIdAndAlgorithms("", list)).collect(Collectors.toList());
        Assert.assertEquals(collapser.collapseAlgorithms(records), expectedResult);
    }

    @Test
    public void testComplexSubtypeAndIntervals() {
        final SVCallRecord cpx1 = new SVCallRecord("cpx1", "chr1", 1000, null,
                "chr1", 1000, null,
                GATKSVVCFConstants.StructuralVariantAnnotationType.CPX,
                GATKSVVCFConstants.ComplexVariantSubtype.dDUP,
                Arrays.asList(new SVCallRecord.ComplexEventInterval("DUP_chr1:6000-8000")),
                null, Collections.singletonList(SVTestUtils.PESR_ALGORITHM),
                Lists.newArrayList(Allele.REF_N, SVTestUtils.CPX_ALLELE),
                Collections.emptyList(), Collections.emptyMap(), Collections.emptySet(), null, SVTestUtils.hg38Dict);
        final SVCallRecord cpx2 = new SVCallRecord("cpx1", "chr1", 1000, null,
                "chr1", 1000, null,
                GATKSVVCFConstants.StructuralVariantAnnotationType.CPX,
                GATKSVVCFConstants.ComplexVariantSubtype.dDUP,
                Arrays.asList(new SVCallRecord.ComplexEventInterval("DUP_chr1:6000-8000")),
                null, Collections.singletonList(SVTestUtils.PESR_ALGORITHM),
                Lists.newArrayList(Allele.REF_N, SVTestUtils.CPX_ALLELE),
                Collections.emptyList(), Collections.emptyMap(), Collections.emptySet(), null, SVTestUtils.hg38Dict);
        final SVCallRecord result = collapser.collapse(new SVClusterEngine.OutputCluster(Lists.newArrayList(cpx1, cpx2)));
        Assert.assertEquals(result.getComplexSubtype(), GATKSVVCFConstants.ComplexVariantSubtype.dDUP);
        Assert.assertEquals(result.getComplexEventIntervals(), cpx1.getComplexEventIntervals());
    }
}