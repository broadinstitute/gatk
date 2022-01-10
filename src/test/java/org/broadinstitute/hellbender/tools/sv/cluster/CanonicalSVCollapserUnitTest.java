package org.broadinstitute.hellbender.tools.sv.cluster;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.StructuralVariantType;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
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

public class CanonicalSVCollapserUnitTest {

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
                        StructuralVariantType.DEL,
                        new Integer[]{0},
                        new Integer[]{0},
                        Collections.emptyList(),
                        Collections.emptyList()
                },
                {
                    Collections.singletonList(Collections.singletonList(Allele.REF_N)),
                        StructuralVariantType.DEL,
                        new Integer[]{1},
                        new Integer[]{1},
                        Collections.emptyList(),
                        Collections.emptyList()
                },
                {
                    Collections.singletonList(Collections.singletonList(Allele.REF_A)),
                        StructuralVariantType.DEL,
                        new Integer[]{1},
                        new Integer[]{1},
                        Collections.emptyList(),
                        Collections.emptyList()
                },
                {
                    Collections.singletonList(Collections.singletonList(Allele.NO_CALL)),
                        StructuralVariantType.DEL,
                        new Integer[]{1},
                        new Integer[]{1},
                        Collections.emptyList(),
                        Collections.emptyList()
                },
                {
                    Collections.singletonList(Collections.singletonList(Allele.ALT_A)),
                        StructuralVariantType.INS,
                        new Integer[]{1},
                        new Integer[]{null},
                        Collections.singletonList(Allele.ALT_A),
                        Collections.singletonList(Allele.ALT_A)
                },
                {
                    Collections.singletonList(Collections.singletonList(Allele.SV_SIMPLE_DEL)),
                        StructuralVariantType.DEL,
                        new Integer[]{1},
                        new Integer[]{0},
                        Collections.singletonList(Allele.SV_SIMPLE_DEL),
                        Collections.singletonList(Allele.SV_SIMPLE_DEL)
                },
                {
                    Lists.newArrayList(Collections.singletonList(Allele.REF_N), Collections.singletonList(Allele.SV_SIMPLE_DEL)),
                        StructuralVariantType.DEL,
                        new Integer[]{1, 1},
                        new Integer[]{1, 0},
                        Collections.singletonList(Allele.SV_SIMPLE_DEL),
                        Collections.singletonList(Allele.SV_SIMPLE_DEL)
                },
                {
                    Lists.newArrayList(Lists.newArrayList(Allele.REF_N, Allele.REF_N), Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL)),
                        StructuralVariantType.DEL,
                        new Integer[]{2, 2},
                        new Integer[]{2, 1},
                        Collections.singletonList(Allele.SV_SIMPLE_DEL),
                        Collections.singletonList(Allele.SV_SIMPLE_DEL)
                },
                {
                    Lists.newArrayList(Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DUP), Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL)),
                        StructuralVariantType.CNV,
                        new Integer[]{2, 2},
                        new Integer[]{3, 1},
                        Lists.newArrayList(Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DUP),
                        Lists.newArrayList(Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DUP)},
                {
                    Lists.newArrayList(Collections.singletonList(MEI_INSERTION_ALLELE), Collections.singletonList(MEI_INSERTION_ALLELE)),
                        StructuralVariantType.INS,
                        new Integer[]{1, 1},
                        new Integer[]{null, null},
                        Collections.singletonList(MEI_INSERTION_ALLELE),
                        Collections.singletonList(MEI_INSERTION_ALLELE)},
                {
                    Lists.newArrayList(Collections.singletonList(Allele.SV_SIMPLE_INS), Collections.singletonList(MEI_INSERTION_ALLELE)),
                        StructuralVariantType.INS,
                        new Integer[]{1, 1},
                        new Integer[]{null, null},
                        Collections.singletonList(Allele.SV_SIMPLE_INS),
                        Collections.singletonList(MEI_INSERTION_ALLELE)},
                {
                    Lists.newArrayList(Collections.singletonList(MEI_INSERTION_ALLELE), Collections.singletonList(SVA_INSERTION_ALLELE)),
                        StructuralVariantType.INS,
                        new Integer[]{1, 1},
                        new Integer[]{null, null},
                        Collections.singletonList(MEI_INSERTION_ALLELE),
                        Collections.singletonList(SVA_INSERTION_ALLELE)},
                {
                    Lists.newArrayList(Collections.singletonList(LINE_INSERTION_ALLELE), Collections.singletonList(SVA_INSERTION_ALLELE)),
                        StructuralVariantType.INS,
                        new Integer[]{1, 1},
                        new Integer[]{null, null},
                        Collections.singletonList(MEI_INSERTION_ALLELE),
                        Collections.singletonList(MEI_INSERTION_ALLELE)},
        };
    }

    @Test(dataProvider= "collapseAltAllelesTestData")
    public void collapseAltAllelesTest(final List<List<Allele>> recordGenotypeAlleles,
                                       final StructuralVariantType svtype,
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

    private static final String TEST_KEY_1 = "TEST_KEY_1";
    private static final String TEST_KEY_2 = "TEST_KEY_2";

    @DataProvider(name = "collapseAttributesTestData")
    public Object[][] collapseAttributesTestData() {
        return new Object[][]{
                // Null value
                {
                        Collections.singletonList("var1"),
                        Collections.singletonList(new String[]{TEST_KEY_1}),
                        Collections.singletonList(new Object[]{null}),
                        2,
                        new String[]{TEST_KEY_1},
                        new Object[]{null}
                },
                // Single key / value
                {
                        Collections.singletonList("var1"),
                        Collections.singletonList(new String[]{TEST_KEY_1}),
                        Collections.singletonList(new Object[]{30}),
                        2,
                        new String[]{TEST_KEY_1},
                        new Object[]{30}
                },
                // Two samples, null values
                {
                        Lists.newArrayList("var1", "var2"),
                        Lists.newArrayList(
                                new String[]{TEST_KEY_1},
                                new String[]{TEST_KEY_1}
                        ),
                        Lists.newArrayList(
                                new Object[]{null},
                                new Object[]{null}),
                        2,
                        new String[]{TEST_KEY_1},
                        new Object[]{null}
                },
                // Two samples, same key/value
                {
                        Lists.newArrayList("var1", "var2"),
                        Lists.newArrayList(
                                new String[]{TEST_KEY_1},
                                new String[]{TEST_KEY_1}
                                ),
                        Lists.newArrayList(
                                new Object[]{30},
                                new Object[]{30}),
                        2,
                        new String[]{TEST_KEY_1},
                        new Object[]{30}
                },
                // Two samples, same key / different value
                {
                        Lists.newArrayList("var1", "var2"),
                        Lists.newArrayList(
                                new String[]{TEST_KEY_1},
                                new String[]{TEST_KEY_1}
                        ),
                        Lists.newArrayList(
                                new Object[]{30},
                                new Object[]{45}),
                        2,
                        new String[]{TEST_KEY_1},
                        new Object[]{null}
                },
                // Two samples, one with an extra key
                {
                        Lists.newArrayList("var1", "var2"),
                        Lists.newArrayList(
                                new String[]{TEST_KEY_1, TEST_KEY_2},
                                new String[]{TEST_KEY_1}
                        ),
                        Lists.newArrayList(
                                new Object[]{30, "VALUE2"},
                                new Object[]{30}),
                        2,
                        new String[]{TEST_KEY_1, TEST_KEY_2},
                        new Object[]{30, "VALUE2"}
                },
        };
    }

    @Test(dataProvider= "collapseAttributesTestData")
    public void collapseAttributesTest(final List<String> variantIds,
                                       final List<String[]> keys,
                                       final List<Object[]> values,
                                       final int expectedCopyNumber,
                                       final String[] expectedKeys,
                                       final Object[] expectedValues) {
        final List<Map<String, Object>> inputAttributesList = IntStream.range(0, keys.size())
                .mapToObj(i -> SVTestUtils.keyValueArraysToMap(keys.get(i), values.get(i)))
                .collect(Collectors.toList());
        final Map<String, Object> expectedAttributes = SVTestUtils.keyValueArraysToMap(expectedKeys, expectedValues);
        expectedAttributes.put(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, expectedCopyNumber);

        // Test as genotype attributes
        final List<Genotype> genotypes = inputAttributesList.stream()
                .map(m -> new GenotypeBuilder().attributes(m).make())
                .collect(Collectors.toList());
        final Map<String, Object> testGenotypeMap = collapser.collapseGenotypeAttributes(genotypes, expectedCopyNumber);
        Assert.assertEquals(testGenotypeMap, expectedAttributes);

        // Test as variant attributes
        final List<SVCallRecord> variants = IntStream.range(0, inputAttributesList.size())
                .mapToObj(i -> SVTestUtils.newNamedDeletionRecordWithAttributes(variantIds.get(i), inputAttributesList.get(i)))
                .collect(Collectors.toList());
        final Map<String, Object> expectedAttributesWithMembers = new HashMap<>(expectedAttributes);
        expectedAttributesWithMembers.put(GATKSVVCFConstants.CLUSTER_MEMBER_IDS_KEY, variantIds);
        expectedAttributesWithMembers.remove(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT);
        final Map<String, Object> testVariantMap = collapser.collapseVariantAttributes(variants);
        Assert.assertEquals(testVariantMap, expectedAttributesWithMembers);
    }

    private Map<String, Object> createGenotypeTestAttributes(final Integer expectedCopyNumber, final Integer copyNumber, final String testVal) {
        final String[] keys = new String[]{
                GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT,
                GATKSVVCFConstants.COPY_NUMBER_FORMAT,
                TEST_KEY_1
        };
        final Object[] vals = new Object[]{expectedCopyNumber, copyNumber, testVal};
        return SVTestUtils.keyValueArraysToMap(keys, vals);
    }

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

    private Map<String, Object> createGenotypeTestAttributes(final Integer expectedCopyNumber) {
        return Collections.singletonMap(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, expectedCopyNumber);
    }

    @DataProvider(name = "collapseSampleGenotypesTestData")
    public Object[][] collapseSampleGenotypesTestData() {
        return new Object[][]{
                // Empty case
                {
                        "sample1",
                        Collections.singletonList(Collections.emptyList()),
                        Collections.singletonList(createGenotypeTestAttributes(0)),
                        Allele.REF_N,
                        Collections.emptyList(),
                        GenotypeBuilder.create(
                                "sample1",
                                Collections.emptyList(),
                                createGenotypeTestAttributes(0)
                        )
                },
                // Extra attribute
                {
                        "sample1",
                        Collections.singletonList(Collections.emptyList()),
                        Collections.singletonList(createGenotypeTestAttributes(0, "test")),
                        Allele.REF_N,
                        Collections.emptyList(),
                        GenotypeBuilder.create(
                                "sample1",
                                Collections.emptyList(),
                                createGenotypeTestAttributes(0, "test")
                        )
                },
                // Haploid no-call
                {
                        "sample1",
                        Lists.newArrayList(
                                Collections.singletonList(Allele.NO_CALL),
                                Collections.singletonList(Allele.NO_CALL)
                        ),
                        Lists.newArrayList(
                                createGenotypeTestAttributes(1),
                                createGenotypeTestAttributes(1)
                        ),
                        Allele.REF_N,
                        Collections.emptyList(),
                        GenotypeBuilder.create(
                                "sample1",
                                Collections.singletonList(Allele.NO_CALL),
                                createGenotypeTestAttributes(1)
                        )
                },
                // Simple ref haploid
                {
                        "sample1",
                        Lists.newArrayList(
                                Collections.singletonList(Allele.REF_N),
                                Collections.singletonList(Allele.REF_N)
                        ),
                        Lists.newArrayList(
                                createGenotypeTestAttributes(1),
                                createGenotypeTestAttributes(1)
                        ),
                        Allele.REF_N,
                        Collections.emptyList(),
                        GenotypeBuilder.create(
                                "sample1",
                                Collections.singletonList(Allele.REF_N),
                                createGenotypeTestAttributes(1)
                        )
                },
                // Simple ref haploid, different ref allele
                {
                        "sample1",
                        Lists.newArrayList(
                                Collections.singletonList(Allele.REF_T),
                                Collections.singletonList(Allele.REF_T)
                        ),
                        Lists.newArrayList(
                                createGenotypeTestAttributes(1),
                                createGenotypeTestAttributes(1)
                        ),
                        Allele.REF_T,
                        Collections.emptyList(),
                        GenotypeBuilder.create(
                                "sample1",
                                Collections.singletonList(Allele.REF_T),
                                createGenotypeTestAttributes(1)
                        )
                },
                // Simple ref haploid, with no-call
                {
                        "sample1",
                        Lists.newArrayList(
                                Collections.singletonList(Allele.REF_T),
                                Collections.singletonList(Allele.NO_CALL)
                        ),
                        Lists.newArrayList(
                                createGenotypeTestAttributes(1),
                                createGenotypeTestAttributes(1)
                        ),
                        Allele.REF_T,
                        Collections.emptyList(),
                        GenotypeBuilder.create(
                                "sample1",
                                Collections.singletonList(Allele.REF_T),
                                createGenotypeTestAttributes(1)
                        )
                },
                // Simple ref diploid
                {
                        "sample1",
                        Lists.newArrayList(
                                Lists.newArrayList(Allele.REF_N, Allele.REF_N),
                                Lists.newArrayList(Allele.REF_N, Allele.REF_N)
                        ),
                        Lists.newArrayList(
                                createGenotypeTestAttributes(2),
                                createGenotypeTestAttributes(2)
                        ),
                        Allele.REF_N,
                        Collections.emptyList(),
                        GenotypeBuilder.create(
                                "sample1",
                                Lists.newArrayList(Allele.REF_N, Allele.REF_N),
                                createGenotypeTestAttributes(2)
                        )
                },
                // Simple ref diploid, with no-call
                {
                        "sample1",
                        Lists.newArrayList(
                                Lists.newArrayList(Allele.REF_N, Allele.REF_N),
                                Lists.newArrayList(Allele.NO_CALL, Allele.NO_CALL)
                        ),
                        Lists.newArrayList(
                                createGenotypeTestAttributes(2),
                                createGenotypeTestAttributes(2)
                        ),
                        Allele.REF_N,
                        Collections.emptyList(),
                        GenotypeBuilder.create(
                                "sample1",
                                Lists.newArrayList(Allele.REF_N, Allele.REF_N),
                                createGenotypeTestAttributes(2)
                        )
                },
                // Simple INS cases
                {
                        "sample1",
                        Lists.newArrayList(
                                Collections.singletonList(Allele.REF_N),
                                Collections.singletonList(Allele.SV_SIMPLE_INS)
                        ),
                        Lists.newArrayList(
                                createGenotypeTestAttributes(1),
                                createGenotypeTestAttributes(1)
                        ),
                        Allele.REF_N,
                        Collections.singletonList(Allele.SV_SIMPLE_INS),
                        GenotypeBuilder.create(
                                "sample1",
                                Lists.newArrayList(Allele.SV_SIMPLE_INS),
                                createGenotypeTestAttributes(1)
                        )
                },
                {
                        "sample1",
                        Lists.newArrayList(
                                Collections.singletonList(Allele.NO_CALL),
                                Collections.singletonList(Allele.SV_SIMPLE_INS)
                        ),
                        Lists.newArrayList(
                                createGenotypeTestAttributes(1),
                                createGenotypeTestAttributes(1)
                        ),
                        Allele.REF_N,
                        Collections.singletonList(Allele.SV_SIMPLE_INS),
                        GenotypeBuilder.create(
                                "sample1",
                                Lists.newArrayList(Allele.SV_SIMPLE_INS),
                                createGenotypeTestAttributes(1)
                        )
                },
                {
                        "sample1",
                        Lists.newArrayList(
                                Collections.singletonList(Allele.SV_SIMPLE_INS),
                                Collections.singletonList(Allele.SV_SIMPLE_INS)
                        ),
                        Lists.newArrayList(
                                createGenotypeTestAttributes(1),
                                createGenotypeTestAttributes(1)
                        ),
                        Allele.REF_N,
                        Collections.singletonList(Allele.SV_SIMPLE_INS),
                        GenotypeBuilder.create(
                                "sample1",
                                Lists.newArrayList(Allele.SV_SIMPLE_INS),
                                createGenotypeTestAttributes(1)
                        )
                },
                // het preferred over hom ref
                {
                        "sample1",
                        Lists.newArrayList(
                                Lists.newArrayList(Allele.REF_N, Allele.REF_N),
                                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_INS)
                        ),
                        Lists.newArrayList(
                                createGenotypeTestAttributes(2),
                                createGenotypeTestAttributes(2)
                        ),
                        Allele.REF_N,
                        Collections.singletonList(Allele.SV_SIMPLE_INS),
                        GenotypeBuilder.create(
                                "sample1",
                                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_INS),
                                createGenotypeTestAttributes(2)
                        )
                },
                // het preferred over hom var
                {
                        "sample1",
                        Lists.newArrayList(
                                Lists.newArrayList(Allele.SV_SIMPLE_INS, Allele.SV_SIMPLE_INS),
                                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_INS)
                        ),
                        Lists.newArrayList(
                                createGenotypeTestAttributes(2),
                                createGenotypeTestAttributes(2)
                        ),
                        Allele.REF_N,
                        Collections.singletonList(Allele.SV_SIMPLE_INS),
                        GenotypeBuilder.create(
                                "sample1",
                                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_INS),
                                createGenotypeTestAttributes(2)
                        )
                },
                // hom more frequent
                {
                        "sample1",
                        Lists.newArrayList(
                                Lists.newArrayList(Allele.SV_SIMPLE_INS, Allele.SV_SIMPLE_INS),
                                Lists.newArrayList(Allele.SV_SIMPLE_INS, Allele.SV_SIMPLE_INS),
                                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_INS)
                        ),
                        Lists.newArrayList(
                                createGenotypeTestAttributes(2),
                                createGenotypeTestAttributes(2),
                                createGenotypeTestAttributes(2)
                        ),
                        Allele.REF_N,
                        Collections.singletonList(Allele.SV_SIMPLE_INS),
                        GenotypeBuilder.create(
                                "sample1",
                                Lists.newArrayList(Allele.SV_SIMPLE_INS, Allele.SV_SIMPLE_INS),
                                createGenotypeTestAttributes(2)
                        )
                },
                // het preferred over both hom ref always and hom var when tied
                {
                        "sample1",
                        Lists.newArrayList(
                                Lists.newArrayList(Allele.SV_SIMPLE_INS, Allele.SV_SIMPLE_INS),
                                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_INS),
                                Lists.newArrayList(Allele.REF_N, Allele.REF_N),
                                Lists.newArrayList(Allele.REF_N, Allele.REF_N)
                        ),
                        Lists.newArrayList(
                                createGenotypeTestAttributes(2),
                                createGenotypeTestAttributes(2),
                                createGenotypeTestAttributes(2),
                                createGenotypeTestAttributes(2)
                        ),
                        Allele.REF_N,
                        Collections.singletonList(Allele.SV_SIMPLE_INS),
                        GenotypeBuilder.create(
                                "sample1",
                                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_INS),
                                createGenotypeTestAttributes(2)
                        )
                },
                // hom is most frequent non-ref
                {
                        "sample1",
                        Lists.newArrayList(
                                Lists.newArrayList(Allele.SV_SIMPLE_INS, Allele.SV_SIMPLE_INS),
                                Lists.newArrayList(Allele.SV_SIMPLE_INS, Allele.SV_SIMPLE_INS),
                                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_INS),
                                Lists.newArrayList(Allele.REF_N, Allele.REF_N),
                                Lists.newArrayList(Allele.REF_N, Allele.REF_N),
                                Lists.newArrayList(Allele.REF_N, Allele.REF_N)
                        ),
                        Lists.newArrayList(
                                createGenotypeTestAttributes(2),
                                createGenotypeTestAttributes(2),
                                createGenotypeTestAttributes(2),
                                createGenotypeTestAttributes(2),
                                createGenotypeTestAttributes(2),
                                createGenotypeTestAttributes(2)
                        ),
                        Allele.REF_N,
                        Collections.singletonList(Allele.SV_SIMPLE_INS),
                        GenotypeBuilder.create(
                                "sample1",
                                Lists.newArrayList(Allele.SV_SIMPLE_INS, Allele.SV_SIMPLE_INS),
                                createGenotypeTestAttributes(2)
                        )
                },
                // triploid - 1 alt
                {
                        "sample1",
                        Lists.newArrayList(
                                Lists.newArrayList(Allele.REF_N, Allele.REF_N, Allele.SV_SIMPLE_INS),
                                Lists.newArrayList(Allele.REF_N, Allele.REF_N, Allele.REF_N)
                        ),
                        Lists.newArrayList(
                                createGenotypeTestAttributes(3),
                                createGenotypeTestAttributes(3)
                        ),
                        Allele.REF_N,
                        Collections.singletonList(Allele.SV_SIMPLE_INS),
                        GenotypeBuilder.create(
                                "sample1",
                                Lists.newArrayList(Allele.REF_N, Allele.REF_N, Allele.SV_SIMPLE_INS),
                                createGenotypeTestAttributes(3)
                        )
                },
                {
                        "sample1",
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
                        Collections.singletonList(Allele.SV_SIMPLE_INS),
                        GenotypeBuilder.create(
                                "sample1",
                                Lists.newArrayList(Allele.REF_N, Allele.REF_N, Allele.SV_SIMPLE_INS),
                                createGenotypeTestAttributes(3)
                        )
                },
                // triploid - 2 alt
                {
                        "sample1",
                        Lists.newArrayList(
                                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_INS, Allele.SV_SIMPLE_INS),
                                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_INS, Allele.SV_SIMPLE_INS),
                                Lists.newArrayList(Allele.REF_N, Allele.REF_N, Allele.SV_SIMPLE_INS),
                                Lists.newArrayList(Allele.REF_N, Allele.REF_N, Allele.REF_N)
                        ),
                        Lists.newArrayList(
                                createGenotypeTestAttributes(3),
                                createGenotypeTestAttributes(3),
                                createGenotypeTestAttributes(3),
                                createGenotypeTestAttributes(3)
                        ),
                        Allele.REF_N,
                        Collections.singletonList(Allele.SV_SIMPLE_INS),
                        GenotypeBuilder.create(
                                "sample1",
                                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_INS, Allele.SV_SIMPLE_INS),
                                createGenotypeTestAttributes(3)
                        )
                },
                // triploid - 3 alt
                {
                        "sample1",
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
                        Collections.singletonList(Allele.SV_SIMPLE_INS),
                        GenotypeBuilder.create(
                                "sample1",
                                Lists.newArrayList(Allele.SV_SIMPLE_INS, Allele.SV_SIMPLE_INS, Allele.SV_SIMPLE_INS),
                                createGenotypeTestAttributes(3)
                        )
                },
                // Simple DEL
                {
                        "sample1",
                        Lists.newArrayList(
                                Collections.singletonList(Allele.REF_N),
                                Collections.singletonList(Allele.SV_SIMPLE_DEL)
                        ),
                        Lists.newArrayList(
                                createGenotypeTestAttributes(1, 1),
                                createGenotypeTestAttributes(1, 0)
                        ),
                        Allele.REF_N,
                        Collections.singletonList(Allele.SV_SIMPLE_DEL),
                        GenotypeBuilder.create(
                                "sample1",
                                Lists.newArrayList(Allele.SV_SIMPLE_DEL),
                                createGenotypeTestAttributes(1, 0)
                        )
                },
                // Simple DEL, falling back on copy number info when GT not available
                {
                        "sample1",
                        Lists.newArrayList(
                                Collections.singletonList(Allele.REF_N),
                                Collections.singletonList(Allele.NO_CALL)
                        ),
                        Lists.newArrayList(
                                createGenotypeTestAttributes(1, 1),
                                createGenotypeTestAttributes(1, 0)
                        ),
                        Allele.REF_N,
                        Collections.singletonList(Allele.SV_SIMPLE_DEL),
                        GenotypeBuilder.create(
                                "sample1",
                                Lists.newArrayList(Allele.SV_SIMPLE_DEL),
                                createGenotypeTestAttributes(1, 0)
                        )
                },
                // Simple DEL diploid het
                {
                        "sample1",
                        Lists.newArrayList(
                                Lists.newArrayList(Allele.REF_N, Allele.REF_N),
                                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL)
                        ),
                        Lists.newArrayList(
                                createGenotypeTestAttributes(2, 2),
                                createGenotypeTestAttributes(2, 1)
                        ),
                        Allele.REF_N,
                        Collections.singletonList(Allele.SV_SIMPLE_DEL),
                        GenotypeBuilder.create(
                                "sample1",
                                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL),
                                createGenotypeTestAttributes(2, 1)
                        )
                },
                // Simple DEL diploid hom var
                {
                        "sample1",
                        Lists.newArrayList(
                                Lists.newArrayList(Allele.REF_N, Allele.REF_N),
                                Lists.newArrayList(Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DEL)
                        ),
                        Lists.newArrayList(
                                createGenotypeTestAttributes(2, 2),
                                createGenotypeTestAttributes(2, 0)
                        ),
                        Allele.REF_N,
                        Collections.singletonList(Allele.SV_SIMPLE_DEL),
                        GenotypeBuilder.create(
                                "sample1",
                                Lists.newArrayList(Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DEL),
                                createGenotypeTestAttributes(2, 0)
                        )
                },
                // Simple DEL diploid hom ref
                {
                        "sample1",
                        Lists.newArrayList(
                                Lists.newArrayList(Allele.REF_N, Allele.REF_N),
                                Lists.newArrayList(Allele.REF_N, Allele.REF_N)
                        ),
                        Lists.newArrayList(
                                createGenotypeTestAttributes(2, 2),
                                createGenotypeTestAttributes(2, 2)
                        ),
                        Allele.REF_N,
                        Collections.singletonList(Allele.SV_SIMPLE_DEL),
                        GenotypeBuilder.create(
                                "sample1",
                                Lists.newArrayList(Allele.REF_N, Allele.REF_N),
                                createGenotypeTestAttributes(2, 2)
                        )
                },
                // Simple DEL triploid with 1 alt
                {
                        "sample1",
                        Lists.newArrayList(
                                Lists.newArrayList(Allele.REF_N, Allele.REF_N, Allele.REF_N),
                                Lists.newArrayList(Allele.REF_N, Allele.REF_N, Allele.SV_SIMPLE_DEL)
                        ),
                        Lists.newArrayList(
                                createGenotypeTestAttributes(3, 3),
                                createGenotypeTestAttributes(3, 2)
                        ),
                        Allele.REF_N,
                        Collections.singletonList(Allele.SV_SIMPLE_DEL),
                        GenotypeBuilder.create(
                                "sample1",
                                Lists.newArrayList(Allele.REF_N, Allele.REF_N, Allele.SV_SIMPLE_DEL),
                                createGenotypeTestAttributes(3, 2)
                        )
                },

                // Simple DUP, haploid
                {
                        "sample1",
                        Lists.newArrayList(
                                Collections.singletonList(Allele.REF_N),
                                Collections.singletonList(Allele.SV_SIMPLE_DUP)
                        ),
                        Lists.newArrayList(
                                createGenotypeTestAttributes(1, 1),
                                createGenotypeTestAttributes(1, 2)
                        ),
                        Allele.REF_N,
                        Collections.singletonList(Allele.SV_SIMPLE_DUP),
                        GenotypeBuilder.create(
                                "sample1",
                                Lists.newArrayList(Allele.SV_SIMPLE_DUP),
                                createGenotypeTestAttributes(1, 2)
                        )
                },
                // Simple DUP, haploid, falling back on copy number info when GT not available and inferring the GT
                {
                        "sample1",
                        Lists.newArrayList(
                                Collections.singletonList(Allele.REF_N),
                                Collections.singletonList(Allele.NO_CALL)
                        ),
                        Lists.newArrayList(
                                createGenotypeTestAttributes(1, 1),
                                createGenotypeTestAttributes(1, 2)
                        ),
                        Allele.REF_N,
                        Collections.singletonList(Allele.SV_SIMPLE_DUP),
                        GenotypeBuilder.create(
                                "sample1",
                                Lists.newArrayList(Allele.SV_SIMPLE_DUP),
                                createGenotypeTestAttributes(1, 2)
                        )
                },
                // Simple DUP, haploid, copy number 3 but phasing is still unambiguous
                {
                        "sample1",
                        Lists.newArrayList(
                                Collections.singletonList(Allele.REF_N),
                                Collections.singletonList(Allele.SV_SIMPLE_DUP)
                        ),
                        Lists.newArrayList(
                                createGenotypeTestAttributes(1, 1),
                                createGenotypeTestAttributes(1, 3)
                        ),
                        Allele.REF_N,
                        Collections.singletonList(Allele.SV_SIMPLE_DUP),
                        GenotypeBuilder.create(
                                "sample1",
                                Lists.newArrayList(Allele.SV_SIMPLE_DUP),
                                createGenotypeTestAttributes(1, 3)
                        )
                },
                // Simple DUP diploid het
                {
                        "sample1",
                        Lists.newArrayList(
                                Lists.newArrayList(Allele.REF_N, Allele.REF_N),
                                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DUP)
                        ),
                        Lists.newArrayList(
                                createGenotypeTestAttributes(2, 2),
                                createGenotypeTestAttributes(2, 3)
                        ),
                        Allele.REF_N,
                        Collections.singletonList(Allele.SV_SIMPLE_DUP),
                        GenotypeBuilder.create(
                                "sample1",
                                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DUP),
                                createGenotypeTestAttributes(2, 3)
                        )
                },
                // Simple DUP diploid hom var - has ambiguous alleles
                {
                        "sample1",
                        Lists.newArrayList(
                                Lists.newArrayList(Allele.REF_N, Allele.REF_N),
                                Lists.newArrayList(Allele.NO_CALL, Allele.NO_CALL)
                        ),
                        Lists.newArrayList(
                                createGenotypeTestAttributes(2, 2),
                                createGenotypeTestAttributes(2, 4)
                        ),
                        Allele.REF_N,
                        Collections.singletonList(Allele.SV_SIMPLE_DUP),
                        GenotypeBuilder.create(
                                "sample1",
                                Lists.newArrayList(Allele.NO_CALL, Allele.NO_CALL),
                                createGenotypeTestAttributes(2, 4)
                        )
                },
                // Simple DUP diploid hom ref
                {
                        "sample1",
                        Lists.newArrayList(
                                Lists.newArrayList(Allele.REF_N, Allele.REF_N),
                                Lists.newArrayList(Allele.REF_N, Allele.REF_N)
                        ),
                        Lists.newArrayList(
                                createGenotypeTestAttributes(2, 2),
                                createGenotypeTestAttributes(2, 2)
                        ),
                        Allele.REF_N,
                        Collections.singletonList(Allele.SV_SIMPLE_DUP),
                        GenotypeBuilder.create(
                                "sample1",
                                Lists.newArrayList(Allele.REF_N, Allele.REF_N),
                                createGenotypeTestAttributes(2, 2)
                        )
                },
                // Simple DUP triploid with 1 alt - unambiguous alleles
                {
                        "sample1",
                        Lists.newArrayList(
                                Lists.newArrayList(Allele.REF_N, Allele.REF_N, Allele.REF_N),
                                Lists.newArrayList(Allele.REF_N, Allele.REF_N, Allele.SV_SIMPLE_DUP)
                        ),
                        Lists.newArrayList(
                                createGenotypeTestAttributes(3, 3),
                                createGenotypeTestAttributes(3, 4)
                        ),
                        Allele.REF_N,
                        Collections.singletonList(Allele.SV_SIMPLE_DUP),
                        GenotypeBuilder.create(
                                "sample1",
                                Lists.newArrayList(Allele.REF_N, Allele.REF_N, Allele.SV_SIMPLE_DUP),
                                createGenotypeTestAttributes(3, 4)
                        )
                },
                // Simple DUP triploid with 2 alts - ambiguous alleles
                {
                        "sample1",
                        Lists.newArrayList(
                                Lists.newArrayList(Allele.REF_N, Allele.REF_N, Allele.REF_N),
                                Lists.newArrayList(Allele.NO_CALL, Allele.NO_CALL, Allele.NO_CALL)
                        ),
                        Lists.newArrayList(
                                createGenotypeTestAttributes(3, 3),
                                createGenotypeTestAttributes(3, 5)
                        ),
                        Allele.REF_N,
                        Collections.singletonList(Allele.SV_SIMPLE_DUP),
                        GenotypeBuilder.create(
                                "sample1",
                                Lists.newArrayList(Allele.NO_CALL, Allele.NO_CALL, Allele.NO_CALL),
                                createGenotypeTestAttributes(3, 5)
                        )
                },
                // Simple DUP triploid where 1-alt genotype should be prioritized
                {
                        "sample1",
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
                        Collections.singletonList(Allele.SV_SIMPLE_DUP),
                        GenotypeBuilder.create(
                                "sample1",
                                Lists.newArrayList(Allele.REF_N, Allele.REF_N, Allele.SV_SIMPLE_DUP),
                                createGenotypeTestAttributes(3, 4)
                        )
                },
                // Simple DUP triploid where 2-alt genotype should be prioritized
                {
                        "sample1",
                        Lists.newArrayList(
                                Lists.newArrayList(Allele.REF_N, Allele.REF_N, Allele.REF_N),
                                Lists.newArrayList(Allele.REF_N, Allele.REF_N, Allele.SV_SIMPLE_DUP),
                                Lists.newArrayList(Allele.NO_CALL, Allele.NO_CALL, Allele.NO_CALL),
                                Lists.newArrayList(Allele.NO_CALL, Allele.NO_CALL, Allele.NO_CALL)
                        ),
                        Lists.newArrayList(
                                createGenotypeTestAttributes(3, 3),
                                createGenotypeTestAttributes(3, 4),
                                createGenotypeTestAttributes(3, 5),
                                createGenotypeTestAttributes(3, 5)
                        ),
                        Allele.REF_N,
                        Collections.singletonList(Allele.SV_SIMPLE_DUP),
                        GenotypeBuilder.create(
                                "sample1",
                                Lists.newArrayList(Allele.NO_CALL, Allele.NO_CALL, Allele.NO_CALL),
                                createGenotypeTestAttributes(3, 5)
                        )
                },


                // Multi-allelic CNV, haploid hom ref
                {
                        "sample1",
                        Lists.newArrayList(
                                Collections.singletonList(Allele.REF_N),
                                Collections.singletonList(Allele.REF_N)
                        ),
                        Lists.newArrayList(
                                createGenotypeTestAttributes(1, 1),
                                createGenotypeTestAttributes(1, 1)
                        ),
                        Allele.REF_N,
                        Lists.newArrayList(Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DUP),
                        GenotypeBuilder.create(
                                "sample1",
                                Lists.newArrayList(Allele.REF_N),
                                createGenotypeTestAttributes(1, 1)
                        )
                },
                // Multi-allelic CNV, haploid del
                {
                        "sample1",
                        Lists.newArrayList(
                                Collections.singletonList(Allele.REF_N),
                                Collections.singletonList(Allele.SV_SIMPLE_DEL)
                        ),
                        Lists.newArrayList(
                                createGenotypeTestAttributes(1, 1),
                                createGenotypeTestAttributes(1, 0)
                        ),
                        Allele.REF_N,
                        Lists.newArrayList(Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DUP),
                        GenotypeBuilder.create(
                                "sample1",
                                Lists.newArrayList(Allele.SV_SIMPLE_DEL),
                                createGenotypeTestAttributes(1, 0)
                        )
                },
                // Multi-allelic CNV, haploid dup
                {
                        "sample1",
                        Lists.newArrayList(
                                Collections.singletonList(Allele.REF_N),
                                Collections.singletonList(Allele.SV_SIMPLE_DUP)
                        ),
                        Lists.newArrayList(
                                createGenotypeTestAttributes(1, 1),
                                createGenotypeTestAttributes(1, 2)
                        ),
                        Allele.REF_N,
                        Lists.newArrayList(Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DUP),
                        GenotypeBuilder.create(
                                "sample1",
                                Lists.newArrayList(Allele.SV_SIMPLE_DUP),
                                createGenotypeTestAttributes(1, 2)
                        )
                },
                // Multi-allelic CNV, diploid hom ref (ambiguous alleles)
                {
                        "sample1",
                        Lists.newArrayList(
                                Lists.newArrayList(Allele.NO_CALL, Allele.NO_CALL),
                                Lists.newArrayList(Allele.NO_CALL, Allele.NO_CALL)
                        ),
                        Lists.newArrayList(
                                createGenotypeTestAttributes(2, 2),
                                createGenotypeTestAttributes(2, 2)
                        ),
                        Allele.REF_N,
                        Lists.newArrayList(Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DUP),
                        GenotypeBuilder.create(
                                "sample1",
                                Lists.newArrayList(Allele.NO_CALL, Allele.NO_CALL),
                                createGenotypeTestAttributes(2, 2)
                        )
                },
                // Multi-allelic CNV, diploid del het
                {
                        "sample1",
                        Lists.newArrayList(
                                Lists.newArrayList(Allele.REF_N, Allele.REF_N),
                                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL)
                        ),
                        Lists.newArrayList(
                                createGenotypeTestAttributes(2, 2),
                                createGenotypeTestAttributes(2, 1)
                        ),
                        Allele.REF_N,
                        Lists.newArrayList(Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DUP),
                        GenotypeBuilder.create(
                                "sample1",
                                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL),
                                createGenotypeTestAttributes(2, 1)
                        )
                },
                // Multi-allelic CNV, diploid del hom
                {
                        "sample1",
                        Lists.newArrayList(
                                Lists.newArrayList(Allele.REF_N, Allele.REF_N),
                                Lists.newArrayList(Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DEL)
                        ),
                        Lists.newArrayList(
                                createGenotypeTestAttributes(2, 2),
                                createGenotypeTestAttributes(2, 0)
                        ),
                        Allele.REF_N,
                        Lists.newArrayList(Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DUP),
                        GenotypeBuilder.create(
                                "sample1",
                                Lists.newArrayList(Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DEL),
                                createGenotypeTestAttributes(2, 0)
                        )
                },
                // Multi-allelic CNV, diploid dup het (ambiguous alleles)
                {
                        "sample1",
                        Lists.newArrayList(
                                Lists.newArrayList(Allele.REF_N, Allele.REF_N),
                                Lists.newArrayList(Allele.NO_CALL, Allele.NO_CALL)
                        ),
                        Lists.newArrayList(
                                createGenotypeTestAttributes(2, 2),
                                createGenotypeTestAttributes(2, 3)
                        ),
                        Allele.REF_N,
                        Lists.newArrayList(Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DUP),
                        GenotypeBuilder.create(
                                "sample1",
                                Lists.newArrayList(Allele.NO_CALL, Allele.NO_CALL),
                                createGenotypeTestAttributes(2, 3)
                        )
                },
                // Multi-allelic CNV, diploid dup hom (ambiguous alleles)
                {
                        "sample1",
                        Lists.newArrayList(
                                Lists.newArrayList(Allele.REF_N, Allele.REF_N),
                                Lists.newArrayList(Allele.NO_CALL, Allele.NO_CALL)
                        ),
                        Lists.newArrayList(
                                createGenotypeTestAttributes(2, 2),
                                createGenotypeTestAttributes(2, 4)
                        ),
                        Allele.REF_N,
                        Lists.newArrayList(Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DUP),
                        GenotypeBuilder.create(
                                "sample1",
                                Lists.newArrayList(Allele.NO_CALL, Allele.NO_CALL),
                                createGenotypeTestAttributes(2, 4)
                        )
                },
                // Multi-allelic CNV, conflicting del and dup genotypes should result in non-call, haploid
                {
                        "sample1",
                        Lists.newArrayList(
                                Collections.singletonList(Allele.SV_SIMPLE_DUP),
                                Collections.singletonList(Allele.SV_SIMPLE_DEL)
                        ),
                        Lists.newArrayList(
                                createGenotypeTestAttributes(1, 2),
                                createGenotypeTestAttributes(1, 0)
                        ),
                        Allele.REF_N,
                        Lists.newArrayList(Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DUP),
                        GenotypeBuilder.create(
                                "sample1",
                                Lists.newArrayList(Allele.NO_CALL),
                                createGenotypeTestAttributes(1, (Integer) null)
                        )
                },
                // Multi-allelic CNV, conflicting del and dup genotypes should result in non-call, diploid
                {
                        "sample1",
                        Lists.newArrayList(
                                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL),
                                Lists.newArrayList(Allele.NO_CALL, Allele.NO_CALL)
                        ),
                        Lists.newArrayList(
                                createGenotypeTestAttributes(2, 1),
                                createGenotypeTestAttributes(2, 3)
                        ),
                        Allele.REF_N,
                        Lists.newArrayList(Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DUP),
                        GenotypeBuilder.create(
                                "sample1",
                                Lists.newArrayList(Allele.NO_CALL, Allele.NO_CALL),
                                createGenotypeTestAttributes(2, (Integer) null)
                        )
                }
        };
    }

    @Test(dataProvider="collapseSampleGenotypesTestData")
    public void collapseSampleGenotypesTest(final String sampleId,
                                            final List<List<Allele>> alleles,
                                            final List<Map<String, Object>> attributes,
                                            final Allele refAllele,
                                            final List<Allele> altAlleles,
                                            final Genotype expected) {
        final List<Genotype> genotypes = IntStream.range(0, alleles.size())
                .mapToObj(i -> GenotypeBuilder.create(sampleId, alleles.get(i), attributes.get(i)))
                .collect(Collectors.toList());
        final Genotype test = collapser.collapseSampleGenotypes(genotypes, refAllele, altAlleles);
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
                {
                        new Integer[]{1000},
                        new String[]{"chr1"},
                        new StructuralVariantType[]{StructuralVariantType.DEL},
                        StructuralVariantType.DEL,
                        null
                },
                {
                        new Integer[]{1000, 1000},
                        new String[]{"chr1", "chr1"},
                        new StructuralVariantType[]{StructuralVariantType.DUP, StructuralVariantType.DUP},
                        StructuralVariantType.DUP,
                        null
                },
                {
                        new Integer[]{300, 400},
                        new String[]{"chr1", "chr1"},
                        new StructuralVariantType[]{StructuralVariantType.DEL, StructuralVariantType.DUP},
                        StructuralVariantType.CNV,
                        null
                },
                {
                        new Integer[]{300, 400},
                        new String[]{"chr1", "chr1"},
                        new StructuralVariantType[]{StructuralVariantType.INS, StructuralVariantType.INS},
                        StructuralVariantType.INS,
                        300
                },
                {
                        new Integer[]{300, 400, 500},
                        new String[]{"chr1", "chr1", "chr1"},
                        new StructuralVariantType[]{StructuralVariantType.INS, StructuralVariantType.INS, StructuralVariantType.INS},
                        StructuralVariantType.INS,
                        400
                },
                {
                        new Integer[]{null},
                        new String[]{"chr2"},
                        new StructuralVariantType[]{StructuralVariantType.BND},
                        StructuralVariantType.BND,
                        null
                },
                {
                        new Integer[]{null, null},
                        new String[]{"chr2", "chr2"},
                        new StructuralVariantType[]{StructuralVariantType.BND, StructuralVariantType.BND},
                        StructuralVariantType.BND,
                        null
                }
        };
    }

    @Test(dataProvider= "collapseLengthTestData")
    public void collapseLengthTest(final Integer[] lengths, final String[] chrom2, final StructuralVariantType[] svtypes,
                                   final StructuralVariantType newType, final Integer expectedLength) {
        final List<SVCallRecord> records = IntStream.range(0, lengths.length).mapToObj(i -> SVTestUtils.newCallRecordWithLengthAndTypeAndChrom2(lengths[i], svtypes[i], chrom2[i])).collect(Collectors.toList());
        Assert.assertEquals(collapser.collapseLength(records, newType), expectedLength);
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
                // 3 variants, ends overlap starts
                {
                        new int[]{1001, 1011, 1021},
                        new int[]{1031, 1011, 1021},
                        StructuralVariantType.DUP,
                        new int[]{1011, 1021},
                        new int[]{1001, 1031},
                        new int[]{1021, 1021}, // min end before max start
                        new int[]{1011, 1021}
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
                        new int[]{1001, 1001},
                        new int[]{1021, 1021},
                        new int[]{1011, 1011}
                },
                // INS, different start/end
                {
                        new int[]{1001, 1011, 1021},
                        new int[]{1011, 1021, 1031},
                        StructuralVariantType.INS,
                        new int[]{1011, 1011},
                        new int[]{1001, 1001},
                        new int[]{1021, 1021},
                        new int[]{1011, 1011}
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

    @Test(expectedExceptions = { IllegalArgumentException.class })
    public void testUnsupportedCNVAltAlleles() throws IllegalArgumentException {
        final List<Allele> siteAltAlleles = Lists.newArrayList(Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_INS);
        CanonicalSVCollapser.getCNVGenotypeAllelesFromCopyNumber(siteAltAlleles, Allele.REF_A, 2, 1);
    }

    @Test(expectedExceptions = { IllegalArgumentException.class })
    public void testUnsupportedCNVAltAllele() throws IllegalArgumentException {
        final List<Allele> siteAltAlleles = Collections.singletonList(Allele.SV_SIMPLE_INS);
        CanonicalSVCollapser.getCNVGenotypeAllelesFromCopyNumber(siteAltAlleles, Allele.REF_A, 2, 2);
    }
}