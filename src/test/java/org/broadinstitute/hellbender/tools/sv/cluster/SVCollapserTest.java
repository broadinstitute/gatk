package org.broadinstitute.hellbender.tools.sv.cluster;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.StructuralVariantType;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SVCallRecordUtils;
import org.broadinstitute.hellbender.tools.sv.SVTestUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class SVCollapserTest {

    private static final CanonicalSVCollapser collapser = new CanonicalSVCollapser(
            SVTestUtils.hg38Reference,
            CanonicalSVCollapser.AltAlleleSummaryStrategy.COMMON_SUBTYPE,
            CanonicalSVCollapser.BreakpointSummaryStrategy.MEDIAN_START_MEDIAN_END,
            CanonicalSVCollapser.InsertionLengthSummaryStrategy.MEDIAN);
    private static final CanonicalSVCollapser collapserMinMax = new CanonicalSVCollapser(
            SVTestUtils.hg38Reference,
            CanonicalSVCollapser.AltAlleleSummaryStrategy.COMMON_SUBTYPE,
            CanonicalSVCollapser.BreakpointSummaryStrategy.MIN_START_MAX_END,
            CanonicalSVCollapser.InsertionLengthSummaryStrategy.MEDIAN);
    private static final CanonicalSVCollapser collapserMaxMin = new CanonicalSVCollapser(
            SVTestUtils.hg38Reference,
            CanonicalSVCollapser.AltAlleleSummaryStrategy.COMMON_SUBTYPE,
            CanonicalSVCollapser.BreakpointSummaryStrategy.MAX_START_MIN_END,
            CanonicalSVCollapser.InsertionLengthSummaryStrategy.MEDIAN);
    private static final CanonicalSVCollapser collapserMean = new CanonicalSVCollapser(
            SVTestUtils.hg38Reference,
            CanonicalSVCollapser.AltAlleleSummaryStrategy.COMMON_SUBTYPE,
            CanonicalSVCollapser.BreakpointSummaryStrategy.MEAN_START_MEAN_END,
            CanonicalSVCollapser.InsertionLengthSummaryStrategy.MEDIAN);
    private static final CanonicalSVCollapser collapserSpecificAltAllele = new CanonicalSVCollapser(
            SVTestUtils.hg38Reference,
            CanonicalSVCollapser.AltAlleleSummaryStrategy.MOST_SPECIFIC_SUBTYPE,
            CanonicalSVCollapser.BreakpointSummaryStrategy.MEDIAN_START_MEDIAN_END,
            CanonicalSVCollapser.InsertionLengthSummaryStrategy.MEDIAN);

    private static final CanonicalSVCollapser collapserInsertionMean = new CanonicalSVCollapser(
            SVTestUtils.hg38Reference,
            CanonicalSVCollapser.AltAlleleSummaryStrategy.COMMON_SUBTYPE,
            CanonicalSVCollapser.BreakpointSummaryStrategy.MEDIAN_START_MEDIAN_END,
            CanonicalSVCollapser.InsertionLengthSummaryStrategy.MEAN);
    private static final CanonicalSVCollapser collapserInsertionMin = new CanonicalSVCollapser(
            SVTestUtils.hg38Reference,
            CanonicalSVCollapser.AltAlleleSummaryStrategy.COMMON_SUBTYPE,
            CanonicalSVCollapser.BreakpointSummaryStrategy.MEDIAN_START_MEDIAN_END,
            CanonicalSVCollapser.InsertionLengthSummaryStrategy.MIN);
    private static final CanonicalSVCollapser collapserInsertionMax = new CanonicalSVCollapser(
            SVTestUtils.hg38Reference,
            CanonicalSVCollapser.AltAlleleSummaryStrategy.COMMON_SUBTYPE,
            CanonicalSVCollapser.BreakpointSummaryStrategy.MEDIAN_START_MEDIAN_END,
            CanonicalSVCollapser.InsertionLengthSummaryStrategy.MAX);
    private static final CanonicalSVCollapser collapserInsertionUndefined = new CanonicalSVCollapser(
            SVTestUtils.hg38Reference,
            CanonicalSVCollapser.AltAlleleSummaryStrategy.COMMON_SUBTYPE,
            CanonicalSVCollapser.BreakpointSummaryStrategy.MEDIAN_START_MEDIAN_END,
            CanonicalSVCollapser.InsertionLengthSummaryStrategy.UNDEFINED);

    private static final CanonicalSVCollapser.AlleleCollectionCollapserComparator alleleComparator = new CanonicalSVCollapser.AlleleCollectionCollapserComparator();

    private static final Allele MEI_INSERTION_ALLELE = Allele.create("<INS:MEI>");
    private static final Allele SVA_INSERTION_ALLELE = Allele.create("<INS:MEI:SVA>");
    private static final Allele LINE_INSERTION_ALLELE = Allele.create("<INS:MEI:LINE>");

    @DataProvider(name = "expectedCopyNumberTestData")
    public Object[][] expectedCopyNumberTestData() {
        return new Object[][]{
                // Result should be same value
                {new int[]{0}, 0},
                {new int[]{1, 1, 1}, 1},
                {new int[]{2, 2, 2}, 2},
                {new int[]{3, 3, 3, 3}, 3}
        };
    }

    @Test(dataProvider= "expectedCopyNumberTestData")
    public void testCollapseExpectedCopyNumber(final int[] input, final int result) {
        final Collection<Genotype> genotypes = IntStream.of(input).mapToObj(p -> SVTestUtils.buildHomGenotypeWithPloidy(Allele.REF_N, p)).collect(Collectors.toList());
        Assert.assertEquals(collapser.collapseExpectedCopyNumber(genotypes), result);
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

    @DataProvider(name = "getSymbolicAlleleSymbolsTestData")
    public Object[][] getSymbolicAlleleSymbolsTestData() {
        return new Object[][]{
                {Allele.create("<DUP>"), new String[] {"DUP"}},
                {Allele.create("<DUP:TANDEM>"), new String[] {"DUP", "TANDEM"}},
                {Allele.create("<INS:MEI:LINE>"), new String[] {"INS", "MEI", "LINE"}}
        };
    }

    @Test(dataProvider= "getSymbolicAlleleSymbolsTestData")
    public void getSymbolicAlleleSymbolsTest(final Allele allele, final String[] result) {
        final String[] test = CanonicalSVCollapser.getSymbolicAlleleSymbols(allele);
        Assert.assertEquals(test, result);
    }

    @DataProvider(name = "collapseAltAllelesTestData")
    public Object[][] collapseAltAllelesTestData() {
        return new Object[][]{
                {Collections.singletonList(Collections.emptyList()), StructuralVariantType.DEL, Collections.emptyList(), Collections.emptyList()},
                {Collections.singletonList(Collections.singletonList(Allele.REF_N)), StructuralVariantType.DEL, Collections.emptyList(), Collections.emptyList()},
                {Collections.singletonList(Collections.singletonList(Allele.REF_A)), StructuralVariantType.DEL, Collections.emptyList(), Collections.emptyList()},
                {Collections.singletonList(Collections.singletonList(Allele.NO_CALL)), StructuralVariantType.DEL, Collections.emptyList(), Collections.emptyList()},
                {Collections.singletonList(Collections.singletonList(Allele.ALT_A)), StructuralVariantType.INS, Collections.singletonList(Allele.ALT_A), Collections.singletonList(Allele.ALT_A)},
                {Collections.singletonList(Collections.singletonList(Allele.SV_SIMPLE_DEL)), StructuralVariantType.DEL, Collections.singletonList(Allele.SV_SIMPLE_DEL), Collections.singletonList(Allele.SV_SIMPLE_DEL)},
                {Lists.newArrayList(Collections.singletonList(Allele.REF_N), Collections.singletonList(Allele.SV_SIMPLE_DEL)), StructuralVariantType.DEL, Collections.singletonList(Allele.SV_SIMPLE_DEL), Collections.singletonList(Allele.SV_SIMPLE_DEL)},
                {Lists.newArrayList(Lists.newArrayList(Allele.REF_N, Allele.REF_N), Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL)), StructuralVariantType.DEL, Collections.singletonList(Allele.SV_SIMPLE_DEL), Collections.singletonList(Allele.SV_SIMPLE_DEL)},
                {Lists.newArrayList(Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DUP), Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL)), StructuralVariantType.CNV, Lists.newArrayList(Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DUP), Lists.newArrayList(Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DUP)},
                {Lists.newArrayList(Collections.singletonList(MEI_INSERTION_ALLELE), Collections.singletonList(MEI_INSERTION_ALLELE)), StructuralVariantType.INS, Collections.singletonList(MEI_INSERTION_ALLELE), Collections.singletonList(MEI_INSERTION_ALLELE)},
                {Lists.newArrayList(Collections.singletonList(Allele.SV_SIMPLE_INS), Collections.singletonList(MEI_INSERTION_ALLELE)), StructuralVariantType.INS, Collections.singletonList(Allele.SV_SIMPLE_INS), Collections.singletonList(MEI_INSERTION_ALLELE)},
                {Lists.newArrayList(Collections.singletonList(MEI_INSERTION_ALLELE), Collections.singletonList(SVA_INSERTION_ALLELE)), StructuralVariantType.INS, Collections.singletonList(MEI_INSERTION_ALLELE), Collections.singletonList(SVA_INSERTION_ALLELE)},
                {Lists.newArrayList(Collections.singletonList(LINE_INSERTION_ALLELE), Collections.singletonList(SVA_INSERTION_ALLELE)), StructuralVariantType.INS, Collections.singletonList(MEI_INSERTION_ALLELE), Collections.singletonList(MEI_INSERTION_ALLELE)},
        };
    }

    @Test(dataProvider= "collapseAltAllelesTestData")
    public void collapseAltAllelesTest(final List<List<Allele>> recordGenotypeAlleles,
                                       final StructuralVariantType svtype,
                                       final List<Allele> resultCommon,
                                       final List<Allele> resultSpecific) {
        final List<Allele> variantAlleles = recordGenotypeAlleles.stream().flatMap(List::stream).distinct().collect(Collectors.toList());
        final List<SVCallRecord> records = recordGenotypeAlleles.stream()
                .map(a -> SVTestUtils.newCallRecordWithAlleles(a, variantAlleles, svtype))
                .collect(Collectors.toList());

        final List<Allele> sortedTestCommon = SVCallRecordUtils.sortAlleles(collapser.collapseAltAlleles(records, svtype));
        final List<Allele> sortedExpectedCommon = SVCallRecordUtils.sortAlleles(resultCommon);
        Assert.assertEquals(sortedTestCommon, sortedExpectedCommon);

        final List<Allele> sortedTestSpecific = SVCallRecordUtils.sortAlleles(collapserSpecificAltAllele.collapseAltAlleles(records, svtype));
        final List<Allele> sortedExpectedSpecific = SVCallRecordUtils.sortAlleles(resultSpecific);
        Assert.assertEquals(sortedTestSpecific, sortedExpectedSpecific);
    }

    @DataProvider(name = "collapseSampleAllelesTestData")
    public Object[][] collapseSampleAllelesTestData() {
        return new Object[][]{
                // empty
                {
                        Collections.singletonList(
                                Collections.emptyList()
                        ),
                        Collections.emptyList(),
                        0,
                        Collections.emptyList()
                },
                // null
                {
                        Collections.singletonList(
                                Collections.singletonList(null)
                        ),
                        Collections.emptyList(),
                        1,
                        Collections.singletonList(Allele.REF_N)
                },
                // REF
                {
                        Collections.singletonList(
                                Collections.singletonList(Allele.REF_N)
                        ),
                        Collections.emptyList(),
                        1,
                        Collections.singletonList(Allele.REF_N)
                },
                // DEL
                {
                    Collections.singletonList(
                            Collections.singletonList(Allele.SV_SIMPLE_DEL)
                    ),
                        Collections.singletonList(Allele.SV_SIMPLE_DEL),
                        1,
                        Collections.singletonList(Allele.SV_SIMPLE_DEL)
                },
                // DEL, DEL
                {
                        Lists.newArrayList(
                                Collections.singletonList(Allele.SV_SIMPLE_DEL),
                                Collections.singletonList(Allele.SV_SIMPLE_DEL)
                        ),
                        Collections.singletonList(Allele.SV_SIMPLE_DEL),
                        1,
                        Collections.singletonList(Allele.SV_SIMPLE_DEL)
                },
                // DEL, REF
                {
                        Lists.newArrayList(
                                Collections.singletonList(Allele.SV_SIMPLE_DEL),
                                Collections.singletonList(Allele.REF_N)
                        ),
                        Collections.singletonList(Allele.SV_SIMPLE_DEL),
                        1,
                        Collections.singletonList(Allele.SV_SIMPLE_DEL)
                },
                // REF/DEL, REF/REF
                {
                        Lists.newArrayList(
                                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL),
                                Lists.newArrayList(Allele.REF_N, Allele.REF_N)
                        ),
                        Collections.singletonList(Allele.SV_SIMPLE_DEL),
                        2,
                        Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL)
                },
                // DEL/DEL, REF/REF
                {
                        Lists.newArrayList(
                                Lists.newArrayList(Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DEL),
                                Lists.newArrayList(Allele.REF_N, Allele.REF_N)
                        ),
                        Collections.singletonList(Allele.SV_SIMPLE_DEL),
                        2,
                        Lists.newArrayList(Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DEL)
                },
                // REF/DEL, REF/DEL
                {
                        Lists.newArrayList(
                                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL),
                                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL)
                        ),
                        Collections.singletonList(Allele.SV_SIMPLE_DEL),
                        2,
                        Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL)
                },
                // REF/DEL, REF/DEL, DEL/DEL
                {
                        Lists.newArrayList(
                                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL),
                                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL),
                                Lists.newArrayList(Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DEL)
                        ),
                        Collections.singletonList(Allele.SV_SIMPLE_DEL),
                        2,
                        Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL)
                },
                // REF/DEL, DEL/DEL, DEL/DEL
                {
                        Lists.newArrayList(
                                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL),
                                Lists.newArrayList(Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DEL),
                                Lists.newArrayList(Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DEL)
                        ),
                        Collections.singletonList(Allele.SV_SIMPLE_DEL),
                        2,
                        Lists.newArrayList(Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DEL)
                },
                // REF/DEL, REF
                {
                        Lists.newArrayList(
                                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL),
                                Collections.singletonList(Allele.REF_N)
                        ),
                        Collections.singletonList(Allele.SV_SIMPLE_DEL),
                        2,
                        Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL)
                },
                // DUP/DEL, REF/REF
                {
                        Lists.newArrayList(
                                Lists.newArrayList(Allele.SV_SIMPLE_DUP, Allele.SV_SIMPLE_DEL),
                                Lists.newArrayList(Allele.REF_N, Allele.REF_N)
                        ),
                        Lists.newArrayList(Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DUP),
                        2,
                        Lists.newArrayList(Allele.SV_SIMPLE_DUP, Allele.SV_SIMPLE_DEL)
                },
        };
    }

    @Test(dataProvider= "collapseSampleAllelesTestData")
    public void collapseSampleAllelesTest(final List<List<Allele>> alleles,
                                          final List<Allele> sampleAltAlleles,
                                          final int expectedCopyNumber,
                                          final List<Allele> result) {
        final Collection<Genotype> genotypes = alleles.stream().map(a -> new GenotypeBuilder().alleles(a).make()).collect(Collectors.toList());
        final List<Allele> sortedTest = SVCallRecordUtils.sortAlleles(collapser.collapseSampleGenotypeAlleles(genotypes, expectedCopyNumber, Allele.REF_N, sampleAltAlleles));
        final List<Allele> sortedResult = SVCallRecordUtils.sortAlleles(result);
        Assert.assertEquals(sortedTest, sortedResult);
    }

    @DataProvider(name = "collapseAttributesTestData")
    public Object[][] collapseAttributesTestData() {
        return new Object[][]{
                // Null value
                {
                        Collections.singletonList("var1"),
                        Collections.singletonList(new String[]{VCFConstants.GENOTYPE_QUALITY_KEY}),
                        Collections.singletonList(new Object[]{null}),
                        Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL),
                        2,
                        new String[]{VCFConstants.GENOTYPE_QUALITY_KEY},
                        new Object[]{null},
                        false
                },
                // Single key / value
                {
                        Collections.singletonList("var1"),
                        Collections.singletonList(new String[]{VCFConstants.GENOTYPE_QUALITY_KEY}),
                        Collections.singletonList(new Object[]{30}),
                        Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL),
                        2,
                        new String[]{VCFConstants.GENOTYPE_QUALITY_KEY},
                        new Object[]{30},
                        false
                },
                // Two samples, null values
                {
                        Lists.newArrayList("var1", "var2"),
                        Lists.newArrayList(
                                new String[]{VCFConstants.GENOTYPE_QUALITY_KEY},
                                new String[]{VCFConstants.GENOTYPE_QUALITY_KEY}
                        ),
                        Lists.newArrayList(
                                new Object[]{null},
                                new Object[]{null}),
                        Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL),
                        2,
                        new String[]{VCFConstants.GENOTYPE_QUALITY_KEY},
                        new Object[]{null},
                        false
                },
                // Two samples, same key/value
                {
                        Lists.newArrayList("var1", "var2"),
                        Lists.newArrayList(
                                new String[]{VCFConstants.GENOTYPE_QUALITY_KEY},
                                new String[]{VCFConstants.GENOTYPE_QUALITY_KEY}
                                ),
                        Lists.newArrayList(
                                new Object[]{30},
                                new Object[]{30}),
                        Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL),
                        2,
                        new String[]{VCFConstants.GENOTYPE_QUALITY_KEY},
                        new Object[]{30},
                        false
                },
                // Two samples, same key / different value
                {
                        Lists.newArrayList("var1", "var2"),
                        Lists.newArrayList(
                                new String[]{VCFConstants.GENOTYPE_QUALITY_KEY},
                                new String[]{VCFConstants.GENOTYPE_QUALITY_KEY}
                        ),
                        Lists.newArrayList(
                                new Object[]{30},
                                new Object[]{45}),
                        Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL),
                        2,
                        new String[]{VCFConstants.GENOTYPE_QUALITY_KEY},
                        new Object[]{null},
                        false
                },
                // Two samples, one with an extra key
                {
                        Lists.newArrayList("var1", "var2"),
                        Lists.newArrayList(
                                new String[]{VCFConstants.GENOTYPE_QUALITY_KEY, "KEY2"},
                                new String[]{VCFConstants.GENOTYPE_QUALITY_KEY}
                        ),
                        Lists.newArrayList(
                                new Object[]{30, "VALUE2"},
                                new Object[]{30}),
                        Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL),
                        2,
                        new String[]{VCFConstants.GENOTYPE_QUALITY_KEY, "KEY2"},
                        new Object[]{30, "VALUE2"},
                        false
                },
                // CNVs
                {
                        Lists.newArrayList("var1", "var2"),
                        Lists.newArrayList(
                                new String[]{GATKSVVCFConstants.COPY_NUMBER_FORMAT},
                                new String[]{GATKSVVCFConstants.COPY_NUMBER_FORMAT}
                        ),
                        Lists.newArrayList(
                                new Object[]{0},
                                new Object[]{0}),
                        Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL),
                        2,
                        new String[]{GATKSVVCFConstants.COPY_NUMBER_FORMAT},
                        new Object[]{1},
                        true
                },
                {
                        Lists.newArrayList("var1", "var2"),
                        Lists.newArrayList(
                                new String[]{GATKSVVCFConstants.COPY_NUMBER_FORMAT},
                                new String[]{GATKSVVCFConstants.COPY_NUMBER_FORMAT}
                        ),
                        Lists.newArrayList(
                                new Object[]{0},
                                new Object[]{0}),
                        Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DUP),
                        2,
                        new String[]{GATKSVVCFConstants.COPY_NUMBER_FORMAT},
                        new Object[]{3},
                        true
                },
                {
                        Lists.newArrayList("var1", "var2"),
                        Lists.newArrayList(
                                new String[]{GATKSVVCFConstants.COPY_NUMBER_FORMAT},
                                new String[]{GATKSVVCFConstants.COPY_NUMBER_FORMAT}
                        ),
                        Lists.newArrayList(
                                new Object[]{0},
                                new Object[]{0}),
                        Lists.newArrayList(Allele.SV_SIMPLE_DUP, Allele.SV_SIMPLE_DEL),
                        2,
                        new String[]{GATKSVVCFConstants.COPY_NUMBER_FORMAT},
                        new Object[]{2},
                        true
                },
                // Edge case with mixed ploidy
                {
                        Lists.newArrayList("var1", "var2"),
                        Lists.newArrayList(
                                new String[]{GATKSVVCFConstants.COPY_NUMBER_FORMAT},
                                new String[]{GATKSVVCFConstants.COPY_NUMBER_FORMAT}
                        ),
                        Lists.newArrayList(
                                new Object[]{2},
                                new Object[]{1}),
                        Lists.newArrayList(Allele.SV_SIMPLE_DEL),
                        1,
                        new String[]{GATKSVVCFConstants.COPY_NUMBER_FORMAT},
                        new Object[]{0},
                        true
                }
        };
    }

    @Test(dataProvider= "collapseAttributesTestData")
    public void collapseAttributesTest(final List<String> variantIds,
                                       final List<String[]> keys,
                                       final List<Object[]> values,
                                       final List<Allele> collapsedAlleles,
                                       final int expectedCopyNumber,
                                       final String[] expectedKeys,
                                       final Object[] expectedValues,
                                       final boolean isCNVTest) {
        final List<Map<String, Object>> inputAttributesList = IntStream.range(0, keys.size())
                .mapToObj(i -> SVTestUtils.keyValueArraysToMap(keys.get(i), values.get(i)))
                .collect(Collectors.toList());
        final Map<String, Object> expectedAttributes = SVTestUtils.keyValueArraysToMap(expectedKeys, expectedValues);
        expectedAttributes.put(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, expectedCopyNumber);

        // Test as genotype attributes
        final List<Genotype> genotypes = inputAttributesList.stream()
                .map(m -> new GenotypeBuilder().attributes(m).make())
                .collect(Collectors.toList());
        Assert.assertEquals(collapser.collapseGenotypeAttributes(genotypes, collapsedAlleles, expectedCopyNumber), expectedAttributes);

        if (!isCNVTest) {
            // Test as variant attributes
            final List<SVCallRecord> variants = IntStream.range(0, inputAttributesList.size())
                    .mapToObj(i -> SVTestUtils.newNamedDeletionRecordWithAttributes(variantIds.get(i), inputAttributesList.get(i)))
                    .collect(Collectors.toList());
            final Map<String, Object> expectedAttributesWithMembers = new HashMap<>(expectedAttributes);
            expectedAttributesWithMembers.put(GATKSVVCFConstants.CLUSTER_MEMBER_IDS_KEY, variantIds);
            expectedAttributesWithMembers.remove(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT);
            Assert.assertEquals(collapser.collapseVariantAttributes(variants), expectedAttributesWithMembers);
        }
    }

    @DataProvider(name = "collapseLengthTestData")
    public Object[][] collapseLengthTestData() {
        return new Object[][]{
                {
                        new Integer[]{1000},
                        new String[]{"chr1"},
                        new StructuralVariantType[]{StructuralVariantType.DEL},
                        1,
                        1000,
                        StructuralVariantType.DEL,
                        1000
                },
                {
                        new Integer[]{1000, 1000},
                        new String[]{"chr1", "chr1"},
                        new StructuralVariantType[]{StructuralVariantType.DUP, StructuralVariantType.DUP},
                        1,
                        1000,
                        StructuralVariantType.DUP,
                        1000
                },
                {
                        new Integer[]{300, 400},
                        new String[]{"chr1", "chr1"},
                        new StructuralVariantType[]{StructuralVariantType.DEL, StructuralVariantType.DUP},
                        1001,
                        1350,
                        StructuralVariantType.CNV,
                        350
                },
                {
                        new Integer[]{300, 400},
                        new String[]{"chr1", "chr1"},
                        new StructuralVariantType[]{StructuralVariantType.INS, StructuralVariantType.INS},
                        1,
                        1,
                        StructuralVariantType.INS,
                        300
                },
                {
                        new Integer[]{300, 400, 500},
                        new String[]{"chr1", "chr1", "chr1"},
                        new StructuralVariantType[]{StructuralVariantType.INS, StructuralVariantType.INS, StructuralVariantType.INS},
                        1,
                        1,
                        StructuralVariantType.INS,
                        400
                },
                {
                        new Integer[]{null},
                        new String[]{"chr2"},
                        new StructuralVariantType[]{StructuralVariantType.BND},
                        1,
                        1,
                        StructuralVariantType.BND,
                        -1
                }
        };
    }

    @Test(dataProvider= "collapseLengthTestData")
    public void collapseLengthTest(final Integer[] lengths, final String[] chrom2, final StructuralVariantType[] svtypes,
                                   final int newStart, final int newEnd, final StructuralVariantType newType, final int expectedLength) {
        final List<SVCallRecord> records = IntStream.range(0, lengths.length).mapToObj(i -> SVTestUtils.newCallRecordWithLengthAndTypeAndChrom2(lengths[i], svtypes[i], chrom2[i])).collect(Collectors.toList());
        Assert.assertEquals(collapser.collapseLength(records, newStart, newEnd, newType), expectedLength);
    }

    @DataProvider(name = "collapseInsertionLengthTestData")
    public Object[][] collapseInsertionLengthTestData() {
        return new Object[][]{
                {
                        new Integer[]{},
                        null,
                        null,
                        null,
                        null,
                        null
                },
                {
                        new Integer[]{null},
                        null,
                        null,
                        null,
                        null,
                        null
                },
                {
                        new Integer[]{100},
                        100,
                        100,
                        100,
                        100,
                        100
                },
                {
                        new Integer[]{null, 100},
                        100,
                        100,
                        100,
                        100,
                        100
                },
                {
                        new Integer[]{100, null},
                        100,
                        100,
                        100,
                        100,
                        100
                },
                {
                        new Integer[]{200, 100},
                        100,
                        150,
                        100,
                        200,
                        null
                },
                {
                        new Integer[]{200, 100, 400},
                        200,
                        234,
                        100,
                        400,
                        null
                }
        };
    }

    @Test(dataProvider= "collapseInsertionLengthTestData")
    public void collapseInsertionLengthTest(final Integer[] lengths,
                                            final Integer expectedMedian,
                                            final Integer expectedMean,
                                            final Integer expectedMin,
                                            final Integer expectedMax,
                                            final Integer expectedUndefined) {
        final List<SVCallRecord> records = IntStream.range(0, lengths.length).mapToObj(i -> SVTestUtils.newCallRecordWithLengthAndTypeAndChrom2(lengths[i], StructuralVariantType.INS, "chr1")).collect(Collectors.toList());
        Assert.assertEquals(collapser.collapseInsertionLength(records), expectedMedian);
        Assert.assertEquals(collapserInsertionMean.collapseInsertionLength(records), expectedMean);
        Assert.assertEquals(collapserInsertionMin.collapseInsertionLength(records), expectedMin);
        Assert.assertEquals(collapserInsertionMax.collapseInsertionLength(records), expectedMax);
        Assert.assertEquals(collapserInsertionUndefined.collapseInsertionLength(records), expectedUndefined);

    }

    @DataProvider(name = "collapseIdsTestData")
    public Object[][] collapseIdsTestData() {
        return new Object[][]{
                {Collections.singletonList("var1"), "var1"},
                {Lists.newArrayList("var1", "var2"), "var1"},
                {Lists.newArrayList("var2", "var1"), "var1"},
        };
    }

    @Test(dataProvider= "collapseIdsTestData")
    public void collapseIdsTest(final List<String> ids, final String expectedResult) {
        final List<SVCallRecord> records = ids.stream().map(SVTestUtils::newDeletionCallRecordWithId).collect(Collectors.toList());
        final String result = collapser.collapseIds(records);
        Assert.assertEquals(result, expectedResult);
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
                                Collections.singletonList("pesr")
                        ),
                        Sets.newHashSet("pesr1")
                },
                {
                        new String[]{"pesr1", "pesr2"},
                        Lists.newArrayList(
                                Collections.singletonList("pesr"),
                                Collections.singletonList("pesr")
                        ),
                        Sets.newHashSet("pesr1", "pesr2")
                },
                {
                        new String[]{"depth1", "pesr1"},
                        Lists.newArrayList(
                                Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
                                Collections.singletonList("pesr")
                        ),
                        Sets.newHashSet("pesr1"),
                },
                {
                        new String[]{"mixed1", "depth1"},
                        Lists.newArrayList(
                                Lists.newArrayList(GATKSVVCFConstants.DEPTH_ALGORITHM, "pesr"),
                                Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM)
                        ),
                        Sets.newHashSet("mixed1"),
                },
                {
                        new String[]{"mixed1", "pesr1"},
                        Lists.newArrayList(
                                Lists.newArrayList(GATKSVVCFConstants.DEPTH_ALGORITHM, "pesr"),
                                Collections.singletonList("pesr")
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
                        new int[]{1001}, // starts
                        new int[]{1100}, // ends
                        StructuralVariantType.DEL,
                        new int[]{1001, 1100}, // median strategy
                        new int[]{1001, 1100}, // min-max strategy
                        new int[]{1001, 1100}, // max-min strategy
                        new int[]{1001, 1100}  // mean strategy
                },
                // 2 variants
                {
                        new int[]{1001, 1011},
                        new int[]{1100, 1110},
                        StructuralVariantType.DEL,
                        new int[]{1001, 1100},
                        new int[]{1001, 1110},
                        new int[]{1011, 1100},
                        new int[]{1006, 1105}
                },
                // 3 variants
                {
                        new int[]{1001, 1011, 1021},
                        new int[]{1100, 1110, 1121},
                        StructuralVariantType.DUP,
                        new int[]{1011, 1110},
                        new int[]{1001, 1121},
                        new int[]{1021, 1100},
                        new int[]{1011, 1110}
                },
                // BND
                {
                        new int[]{1001, 1011, 1021},
                        new int[]{1100, 1110, 1121},
                        StructuralVariantType.BND,
                        new int[]{1011, 1110},
                        new int[]{1001, 1121},
                        new int[]{1021, 1100},
                        new int[]{1011, 1110}
                },
                // INS, same start/end
                {
                        new int[]{1001, 1011, 1021},
                        new int[]{1001, 1011, 1021},
                        StructuralVariantType.INS,
                        new int[]{1011, 1011},
                        new int[]{1011, 1011},
                        new int[]{1011, 1011},
                        new int[]{1011, 1011}
                },
                // INS, different start/end
                {
                        new int[]{1001, 1011, 1021},
                        new int[]{1011, 1021, 1031},
                        StructuralVariantType.INS,
                        new int[]{1016, 1016},
                        new int[]{1016, 1016},
                        new int[]{1016, 1016},
                        new int[]{1016, 1016}
                }
        };
    }

    @Test(dataProvider= "collapseIntervalTestData")
    public void collapseIntervalTest(final int[] starts, final int[] ends, final StructuralVariantType svtype,
                                     final int[] expectedMedian, final int[] expectedMinMax, final int[] expectedMaxMin,
                                     final int[] expectedMean) {
        final List<SVCallRecord> records =  IntStream.range(0, starts.length)
                .mapToObj(i -> SVTestUtils.newCallRecordWithIntervalAndType(starts[i], ends[i], svtype)).collect(Collectors.toList());

        final Pair<Integer, Integer> resultMedian = collapser.collapseInterval(records);
        Assert.assertEquals((int) resultMedian.getKey(), expectedMedian[0]);
        Assert.assertEquals((int) resultMedian.getValue(), expectedMedian[1]);

        final Pair<Integer, Integer> resultMinMax = collapserMinMax.collapseInterval(records);
        Assert.assertEquals((int) resultMinMax.getKey(), expectedMinMax[0]);
        Assert.assertEquals((int) resultMinMax.getValue(), expectedMinMax[1]);

        final Pair<Integer, Integer> resultMaxMin = collapserMaxMin.collapseInterval(records);
        Assert.assertEquals((int) resultMaxMin.getKey(), expectedMaxMin[0]);
        Assert.assertEquals((int) resultMaxMin.getValue(), expectedMaxMin[1]);

        final Pair<Integer, Integer> resultMean = collapserMean.collapseInterval(records);
        Assert.assertEquals((int) resultMean.getKey(), expectedMean[0]);
        Assert.assertEquals((int) resultMean.getValue(), expectedMean[1]);
    }

    @DataProvider(name = "collapseTypesTestData")
    public Object[][] collapseTypesTestData() {
        return new Object[][]{
                {Collections.singletonList(StructuralVariantType.DEL), StructuralVariantType.DEL},
                {Collections.singletonList(StructuralVariantType.DUP), StructuralVariantType.DUP},
                {Collections.singletonList(StructuralVariantType.INS), StructuralVariantType.INS},
                {Collections.singletonList(StructuralVariantType.INV), StructuralVariantType.INV},
                {Collections.singletonList(StructuralVariantType.BND), StructuralVariantType.BND},
                {Collections.singletonList(StructuralVariantType.CNV), StructuralVariantType.CNV},
                {Lists.newArrayList(StructuralVariantType.DEL, StructuralVariantType.DUP), StructuralVariantType.CNV},
                {Lists.newArrayList(StructuralVariantType.DEL, StructuralVariantType.DEL), StructuralVariantType.DEL},
                {Lists.newArrayList(StructuralVariantType.DUP, StructuralVariantType.DUP), StructuralVariantType.DUP},
                {Lists.newArrayList(StructuralVariantType.INS, StructuralVariantType.INS), StructuralVariantType.INS},
                {Lists.newArrayList(StructuralVariantType.INV, StructuralVariantType.INV), StructuralVariantType.INV},
                {Lists.newArrayList(StructuralVariantType.BND, StructuralVariantType.BND), StructuralVariantType.BND},
                {Lists.newArrayList(StructuralVariantType.CNV, StructuralVariantType.CNV), StructuralVariantType.CNV},
                {Lists.newArrayList(StructuralVariantType.DEL, StructuralVariantType.DUP, StructuralVariantType.CNV), StructuralVariantType.CNV}
        };
    }

    @Test(dataProvider= "collapseTypesTestData")
    public void collapseTypesTest(final List<StructuralVariantType> types, final StructuralVariantType expectedResult) {
        final List<SVCallRecord> records = types.stream().map(t -> SVTestUtils.newCallRecordWithIntervalAndType(1, 100, t)).collect(Collectors.toList());
        Assert.assertEquals(collapser.collapseTypes(records), expectedResult);
    }

    @DataProvider(name = "collapseAlgorithmsTestData")
    public Object[][] collapseAlgorithmsTestData() {
        return new Object[][]{
                {Collections.singletonList(Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM)), Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM)},
                {Collections.singletonList(Collections.singletonList("pesr")), Collections.singletonList("pesr")},
                {Collections.singletonList(Lists.newArrayList("pesrA", "pesrB")), Lists.newArrayList("pesrA", "pesrB")},
                {Collections.singletonList(Lists.newArrayList("pesrB", "pesrA")), Lists.newArrayList("pesrA", "pesrB")},
                {Lists.newArrayList(Collections.singletonList("pesrA"), Collections.singletonList("pesrB")), Lists.newArrayList("pesrA", "pesrB")},
                {Lists.newArrayList(Collections.singletonList("pesrB"), Collections.singletonList("pesrA")), Lists.newArrayList("pesrA", "pesrB")},
                {Lists.newArrayList(Lists.newArrayList("pesrA", "pesrB"), Collections.singletonList("pesrB")), Lists.newArrayList("pesrA", "pesrB")},
                {Lists.newArrayList(Lists.newArrayList("pesrB", "pesrA"), Collections.singletonList("pesrA")), Lists.newArrayList("pesrA", "pesrB")},
        };
    }

    @Test(dataProvider= "collapseAlgorithmsTestData")
    public void collapseAlgorithmsTest(final List<List<String>> algorithmLists, final List<String> expectedResult) {
        final List<SVCallRecord> records = algorithmLists.stream().map(list -> SVTestUtils.newDeletionCallRecordWithIdAndAlgorithms("", list)).collect(Collectors.toList());
        Assert.assertEquals(collapser.collapseAlgorithms(records), expectedResult);
    }

    @DataProvider(name = "alleleCollapserComparatorTestData")
    public Object[][] alleleCollapserComparatorTestData() {
        return new Object[][]{
                {Collections.emptyList(), Collections.emptyList(), 0},
                {Collections.singletonList(Allele.NO_CALL), Collections.singletonList(Allele.NO_CALL), 0},
                {Collections.singletonList(Allele.REF_N), Collections.singletonList(Allele.REF_N), 0},
                {Collections.singletonList(Allele.SV_SIMPLE_DEL), Collections.singletonList(Allele.SV_SIMPLE_DEL), 0},
                {Lists.newArrayList(Allele.SV_SIMPLE_DUP, Allele.SV_SIMPLE_DUP), Lists.newArrayList(Allele.SV_SIMPLE_DUP, Allele.SV_SIMPLE_DUP), 0},
                // When otherwise equal up to common length, shorter list first should return 1 <=> longer list first should return -1
                {Lists.newArrayList(Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DEL), Collections.singletonList(Allele.SV_SIMPLE_DEL), 1},
                {Collections.singletonList(Allele.SV_SIMPLE_DEL), Lists.newArrayList(Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DUP), -1},
                {Lists.newArrayList(Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DUP), Collections.singletonList(Allele.SV_SIMPLE_DEL), 1},
                {Collections.singletonList(Allele.SV_SIMPLE_DEL), Lists.newArrayList(Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DUP), -1},
                {Lists.newArrayList(Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DUP), Lists.newArrayList(Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DEL), 1},
                {Lists.newArrayList(Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DEL), Lists.newArrayList(Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DUP), -1}
        };
    }

    @Test(dataProvider= "alleleCollapserComparatorTestData")
    public void alleleCollapserComparatorTest(final List<Allele> allelesA, final List<Allele> allelesB, final int expectedResult) {
        Assert.assertEquals(alleleComparator.compare(allelesA, allelesB), expectedResult);
    }
}