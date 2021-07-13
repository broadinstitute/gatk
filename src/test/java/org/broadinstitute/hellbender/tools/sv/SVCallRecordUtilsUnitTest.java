package org.broadinstitute.hellbender.tools.sv;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.VCFConstants;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;
import java.util.stream.Collectors;

public class SVCallRecordUtilsUnitTest {

    private static final Allele BND_ALLELE = Allele.create("<BND>");
    private static final List<Allele> ALLELES_DEL = Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL);
    private static final List<Allele> ALLELES_INS = Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_INS);
    private static final List<Allele> ALLELES_BND = Lists.newArrayList(Allele.REF_N, BND_ALLELE);

    private static final List<String> DEPTH_ONLY_ALGORITHM = Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM);
    private static final List<String> PESR_ALGORITHM = Collections.singletonList("pesr");

    private static final Genotype GENOTYPE_DEL_1 = new GenotypeBuilder("sample1")
            .alleles(Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL)).make();
    private static final Genotype GENOTYPE_DEL_2 = new GenotypeBuilder("sample2")
            .alleles(Lists.newArrayList(Allele.REF_N, Allele.REF_N)).make();
    private static final Genotype GENOTYPE_DEL_3 = new GenotypeBuilder("sample3")
            .alleles(Lists.newArrayList(Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DEL)).make();

    private static final Genotype GENOTYPE_INS_1 = new GenotypeBuilder("sample1")
            .alleles(Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_INS)).make();
    private static final Genotype GENOTYPE_INS_2 = new GenotypeBuilder("sample2")
            .alleles(Lists.newArrayList(Allele.SV_SIMPLE_INS, Allele.SV_SIMPLE_INS)).make();

    private static final Genotype GENOTYPE_BND_1 = new GenotypeBuilder("sample1")
            .alleles(Lists.newArrayList(Allele.REF_N, BND_ALLELE)).make();

    private static final Comparator<SVCallRecord> RECORD_COMPARATOR = SVCallRecordUtils.getCallComparator(SVTestUtils.hg38Dict);

    @DataProvider(name = "testGetVariantBuilderData")
    public Object[][] testGetVariantBuilderData() {
        return new Object[][]{
                // DEL
                {
                        new SVCallRecord("var1", "chr1", 1000, true, "chr1", 1999, false, StructuralVariantType.DEL, 1000,
                            DEPTH_ONLY_ALGORITHM,
                            ALLELES_DEL,
                            Lists.newArrayList(GENOTYPE_DEL_1, GENOTYPE_DEL_2),
                            Collections.emptyMap()),
                        new VariantContextBuilder("", "chr1", 1000, 1999, ALLELES_DEL)
                                .id("var1")
                                .genotypes(GENOTYPE_DEL_1, GENOTYPE_DEL_2)
                                .attribute(VCFConstants.END_KEY, 1999)
                                .attribute(GATKSVVCFConstants.ALGORITHMS_ATTRIBUTE, DEPTH_ONLY_ALGORITHM)
                                .attribute(GATKSVVCFConstants.SVLEN, 1000)
                                .attribute(GATKSVVCFConstants.SVTYPE, StructuralVariantType.DEL)
                                .make(),
                        Collections.emptyList()
                },
                // DEL w/ null ref allele
                {
                        new SVCallRecord("var1", "chr1", 1000, true, "chr1", 1999, false, StructuralVariantType.DEL, 1000,
                                DEPTH_ONLY_ALGORITHM,
                                Collections.singletonList(Allele.SV_SIMPLE_DEL),
                                Collections.singletonList(GENOTYPE_DEL_3),
                                Collections.emptyMap()),
                        new VariantContextBuilder("", "chr1", 1000, 1999, ALLELES_DEL)
                                .id("var1")
                                .genotypes(GENOTYPE_DEL_3)
                                .attribute(VCFConstants.END_KEY, 1999)
                                .attribute(GATKSVVCFConstants.ALGORITHMS_ATTRIBUTE, DEPTH_ONLY_ALGORITHM)
                                .attribute(GATKSVVCFConstants.SVLEN, 1000)
                                .attribute(GATKSVVCFConstants.SVTYPE, StructuralVariantType.DEL)
                                .make(),
                        Collections.emptyList()
                },
                // INS
                {
                        new SVCallRecord("var2", "chr1", 1000, true, "chr1", 1000, false, StructuralVariantType.INS, 500,
                                PESR_ALGORITHM,
                                ALLELES_INS,
                                Lists.newArrayList(GENOTYPE_INS_1),
                                Collections.emptyMap()),
                        new VariantContextBuilder("", "chr1", 1000, 1000, ALLELES_INS)
                                .id("var2")
                                .genotypes(GENOTYPE_INS_1)
                                .attribute(VCFConstants.END_KEY, 1000)
                                .attribute(GATKSVVCFConstants.ALGORITHMS_ATTRIBUTE, PESR_ALGORITHM)
                                .attribute(GATKSVVCFConstants.SVLEN, 500)
                                .attribute(GATKSVVCFConstants.SVTYPE, StructuralVariantType.INS)
                                .make(),
                        Collections.emptyList()
                },
                // INS, undefined length
                {
                        new SVCallRecord("var2", "chr1", 1000, true, "chr1", 1000, false, StructuralVariantType.INS, null,
                                PESR_ALGORITHM,
                                ALLELES_INS,
                                Lists.newArrayList(GENOTYPE_INS_1),
                                Collections.emptyMap()),
                        new VariantContextBuilder("", "chr1", 1000, 1000, ALLELES_INS)
                                .id("var2")
                                .genotypes(GENOTYPE_INS_1)
                                .attribute(VCFConstants.END_KEY, 1000)
                                .attribute(GATKSVVCFConstants.ALGORITHMS_ATTRIBUTE, PESR_ALGORITHM)
                                .attribute(GATKSVVCFConstants.SVTYPE, StructuralVariantType.INS)
                                .make(),
                        Collections.emptyList()
                },
                // BND
                {
                        new SVCallRecord("var3", "chr1", 1000, false, "chr2", 1999, true, StructuralVariantType.BND, null,
                                PESR_ALGORITHM,
                                ALLELES_BND,
                                Lists.newArrayList(GENOTYPE_BND_1),
                                Collections.emptyMap()),
                        new VariantContextBuilder("", "chr1", 1000, 1000, ALLELES_BND)
                                .id("var3")
                                .genotypes(GENOTYPE_BND_1)
                                .attribute(VCFConstants.END_KEY, 1000)
                                .attribute(GATKSVVCFConstants.ALGORITHMS_ATTRIBUTE, PESR_ALGORITHM)
                                .attribute(GATKSVVCFConstants.STRANDS_ATTRIBUTE, "-+")
                                .attribute(GATKSVVCFConstants.SVTYPE, StructuralVariantType.BND)
                                .attribute(GATKSVVCFConstants.CONTIG2_ATTRIBUTE, "chr2")
                                .attribute(GATKSVVCFConstants.END2_ATTRIBUTE, 1999)
                                .make(),
                        Collections.emptyList()
                },
        };
    }

    @Test(dataProvider= "testGetVariantBuilderData")
    public void testGetVariantBuilder(final SVCallRecord record, final VariantContext expected, final List<String> attributesToIgnore) {
        final VariantContext result = SVCallRecordUtils.getVariantBuilder(record).make();
        VariantContextTestUtils.assertVariantContextsAreEqual(result, expected, attributesToIgnore, Collections.emptyList());
    }

    @Test
    public void testFillMissingSamplesWithGenotypes() {
        final ArrayList<Genotype> genotypes = new ArrayList<>();
        genotypes.add(GENOTYPE_DEL_1);
        genotypes.add(GENOTYPE_DEL_2);
        final GenotypesContext genotypesContext = GenotypesContext.create(genotypes);
        final Set<String> samplesNoMissing = Sets.newHashSet(GENOTYPE_DEL_1.getSampleName(), GENOTYPE_DEL_2.getSampleName());
        final Set<String> samples = Sets.newHashSet(GENOTYPE_DEL_1.getSampleName(), GENOTYPE_DEL_2.getSampleName(),
                "sample3", "sample4");
        final List<Allele> alleles = Lists.newArrayList(Allele.NO_CALL);
        final Map<String, Object> attributes = new HashMap<>();
        attributes.put(GATKSVVCFConstants.ALGORITHMS_ATTRIBUTE, Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM));

        final GenotypesContext resultNoMissing = SVCallRecordUtils.populateGenotypesForMissingSamplesWithAlleles(genotypesContext, samplesNoMissing, alleles, attributes);
        Assert.assertEquals(resultNoMissing, genotypesContext);

        final GenotypesContext result = SVCallRecordUtils.populateGenotypesForMissingSamplesWithAlleles(genotypesContext, samples, alleles, attributes);
        Assert.assertEquals(result.size(), 4);

        final Genotype g1 = result.get(GENOTYPE_DEL_1.getSampleName());
        final Genotype g2 = result.get(GENOTYPE_DEL_2.getSampleName());
        final Genotype g3 = result.get("sample3");
        final Genotype g4 = result.get("sample4");

        final Genotype g3Expected = new GenotypeBuilder("sample3", alleles).attributes(attributes).make();
        final Genotype g4Expected = new GenotypeBuilder("sample4", alleles).attributes(attributes).make();

        VariantContextTestUtils.assertGenotypesAreEqual(g1, GENOTYPE_DEL_1);
        VariantContextTestUtils.assertGenotypesAreEqual(g2, GENOTYPE_DEL_2);
        VariantContextTestUtils.assertGenotypesAreEqual(g3, g3Expected);
        VariantContextTestUtils.assertGenotypesAreEqual(g4, g4Expected);
    }

    @Test
    public void testCopyCallWithNewGenotypes() {
        final SVCallRecord record = new SVCallRecord("var1", "chr1", 1000, true, "chr1", 1999, false, StructuralVariantType.DEL, 1000,
                DEPTH_ONLY_ALGORITHM,
                ALLELES_DEL,
                Lists.newArrayList(GENOTYPE_DEL_1, GENOTYPE_DEL_2),
                Collections.singletonMap(GATKSVVCFConstants.COPY_NUMBER_FORMAT, 0));
        final GenotypesContext genotypes = GenotypesContext.copy(Collections.singletonList(GENOTYPE_DEL_3));
        final SVCallRecord result = SVCallRecordUtils.copyCallWithNewGenotypes(record, genotypes);
        final SVCallRecord expected = new SVCallRecord("var1", "chr1", 1000, true, "chr1", 1999, false, StructuralVariantType.DEL, 1000,
                DEPTH_ONLY_ALGORITHM,
                ALLELES_DEL,
                genotypes,
                Collections.singletonMap(GATKSVVCFConstants.COPY_NUMBER_FORMAT, 0));
        SVTestUtils.assertEqualsExceptMembership(result, expected);
    }

    @DataProvider(name = "testGetCallComparatorData")
    public Object[][] testGetCallComparatorData() {
        return new Object[][]{
                // Identical records
                {
                        SVTestUtils.newCallRecordWithCoordinates("var1", "chr1", 1000, "chr1", 1999),
                        SVTestUtils.newCallRecordWithCoordinates("var2", "chr1", 1000, "chr1", 1999),
                        0
                },
                // Position
                {
                        SVTestUtils.newCallRecordWithCoordinates("var1", "chr1", 1000, "chr1", 1999),
                        SVTestUtils.newCallRecordWithCoordinates("var2", "chr1", 1001, "chr1", 1999),
                        -1
                },
                {
                        SVTestUtils.newCallRecordWithCoordinates("var1", "chr1", 1001, "chr1", 1999),
                        SVTestUtils.newCallRecordWithCoordinates("var2", "chr1", 1000, "chr1", 1999),
                        1
                },
                {
                        SVTestUtils.newCallRecordWithCoordinates("var1", "chr1", 1000, "chr1", 1999),
                        SVTestUtils.newCallRecordWithCoordinatesAndType("var2", "chrX", 1000, "chr1", 1999, StructuralVariantType.BND),
                        -1
                },
                {
                        SVTestUtils.newCallRecordWithCoordinatesAndType("var1", "chrX", 1000, "chr1", 1999, StructuralVariantType.BND),
                        SVTestUtils.newCallRecordWithCoordinates("var2", "chr1", 1000, "chr1", 1999),
                        1
                },
                {
                        SVTestUtils.newCallRecordWithCoordinates("var1", "chr1", 1000, "chr1", 1999),
                        SVTestUtils.newCallRecordWithCoordinates("var2", "chr1", 1000, "chr1", 2000),
                        -1
                },
                {
                        SVTestUtils.newCallRecordWithCoordinates("var1", "chr1", 1000, "chr1", 2000),
                        SVTestUtils.newCallRecordWithCoordinates("var2", "chr1", 1000, "chr1", 1999),
                        1
                },
                {
                        SVTestUtils.newCallRecordWithCoordinates("var1", "chr1", 1000, "chr1", 1999),
                        SVTestUtils.newCallRecordWithCoordinatesAndType("var2", "chr1", 1000, "chrX", 1999, StructuralVariantType.BND),
                        -1
                },
                {
                        SVTestUtils.newCallRecordWithCoordinatesAndType("var1", "chr1", 1000, "chrX", 1999, StructuralVariantType.BND),
                        SVTestUtils.newCallRecordWithCoordinates("var2", "chr1", 1000, "chr1", 1999),
                        1
                },
                // Strands
                {
                        SVTestUtils.newBndCallRecordWithStrands(true, false),
                        SVTestUtils.newBndCallRecordWithStrands(false, false),
                        Boolean.compare(true, false)
                },
                {
                        SVTestUtils.newBndCallRecordWithStrands(false, false),
                        SVTestUtils.newBndCallRecordWithStrands(true, false),
                        Boolean.compare(false, true)
                },
                {
                        SVTestUtils.newBndCallRecordWithStrands(true, true),
                        SVTestUtils.newBndCallRecordWithStrands(true, false),
                        Boolean.compare(true, false)
                },
                {
                        SVTestUtils.newBndCallRecordWithStrands(true, false),
                        SVTestUtils.newBndCallRecordWithStrands(true, true),
                        Boolean.compare(false, true)
                },
                {
                        SVTestUtils.newCnvCallRecordWithStrands(null, null),
                        SVTestUtils.newCnvCallRecordWithStrands(null, null),
                        0
                },
                // Length
                {
                        SVTestUtils.newCallRecordInsertionWithLength(null),
                        SVTestUtils.newCallRecordInsertionWithLength(null),
                        0
                },
                {
                        SVTestUtils.newCallRecordInsertionWithLength(100),
                        SVTestUtils.newCallRecordInsertionWithLength(200),
                        -1
                },
                {
                        SVTestUtils.newCallRecordInsertionWithLength(200),
                        SVTestUtils.newCallRecordInsertionWithLength(100),
                        1
                },
                {
                        SVTestUtils.newCallRecordInsertionWithLength(null),
                        SVTestUtils.newCallRecordInsertionWithLength(100),
                        -1
                },
                {
                        SVTestUtils.newCallRecordInsertionWithLength(100),
                        SVTestUtils.newCallRecordInsertionWithLength(null),
                        1
                },
                // SV type
                {
                        SVTestUtils.newCallRecordWithIntervalAndType(1000, 1999, StructuralVariantType.DEL),
                        SVTestUtils.newCallRecordWithIntervalAndType(1000, 1999, StructuralVariantType.DUP),
                        StructuralVariantType.DEL.compareTo(StructuralVariantType.DUP)
                },
                {
                        SVTestUtils.newCallRecordWithIntervalAndType(1000, 1999, StructuralVariantType.DUP),
                        SVTestUtils.newCallRecordWithIntervalAndType(1000, 1999, StructuralVariantType.DEL),
                        StructuralVariantType.DUP.compareTo(StructuralVariantType.DEL)
                },
                {
                        SVTestUtils.newCallRecordWithIntervalAndType(1000, 1000, StructuralVariantType.BND),
                        SVTestUtils.newCallRecordWithIntervalAndType(1000, 1000, StructuralVariantType.INS),
                        StructuralVariantType.BND.compareTo(StructuralVariantType.INS)
                }
        };
    }

    @Test(dataProvider= "testGetCallComparatorData")
    public void testGetCallComparator(final SVCallRecord record1, final SVCallRecord record2, final int expected) {
        Assert.assertEquals(RECORD_COMPARATOR.compare(record1, record2), expected);
    }

    @Test
    public void testConvertInversionsToBreakends() {
        final SVCallRecord nonInversion = SVTestUtils.newCallRecordWithIntervalAndType(1000, 1999, StructuralVariantType.DEL);
        final List<SVCallRecord> nonInversionResult = SVCallRecordUtils.convertInversionsToBreakends(nonInversion).collect(Collectors.toList());
        Assert.assertEquals(nonInversionResult.size(), 1);
        Assert.assertNotNull(nonInversionResult.get(0));
        SVTestUtils.assertEqualsExceptMembership(nonInversionResult.get(0), nonInversion);

        final SVCallRecord inversion = new SVCallRecord("", "chr1", 1000, true, "chr1", 1999, true, StructuralVariantType.INV, 1000,
                Collections.singletonList("pesr"),
                Collections.emptyList(),
                Collections.emptyList(),
                Collections.emptyMap());
        final List<SVCallRecord> inversionResult = SVCallRecordUtils.convertInversionsToBreakends(inversion).collect(Collectors.toList());
        Assert.assertEquals(inversionResult.size(), 2);

        final SVCallRecord bnd1 = inversionResult.get(0);
        final SVCallRecord bnd2 = inversionResult.get(1);
        Assert.assertNotNull(bnd1);
        Assert.assertNotNull(bnd2);
        Assert.assertEquals(bnd1.getType(), StructuralVariantType.BND);
        Assert.assertEquals(bnd2.getType(), StructuralVariantType.BND);
        Assert.assertNull(bnd1.getLength());
        Assert.assertNull(bnd2.getLength());
        Assert.assertEquals(bnd1.getStrandA(), bnd1.getStrandB());
        Assert.assertEquals(bnd2.getStrandA(), bnd2.getStrandB());
        Assert.assertNotEquals(bnd1.getStrandA(), bnd2.getStrandA());
        Assert.assertEquals(bnd1.getPositionA(), 1000);
        Assert.assertEquals(bnd1.getPositionB(), 1999);
        Assert.assertEquals(bnd2.getPositionA(), 1000);
        Assert.assertEquals(bnd2.getPositionB(), 1999);
    }

    @DataProvider(name = "testIsAltGenotypeData")
    public Object[][] testIsAltGenotypeData() {
        return new Object[][]{
                {new GenotypeBuilder(), false},
                {new GenotypeBuilder().alleles(Collections.singletonList(null)), false},
                {new GenotypeBuilder().alleles(Collections.singletonList(Allele.REF_N)), false},
                {new GenotypeBuilder().alleles(Collections.singletonList(Allele.NO_CALL)), false},
                {new GenotypeBuilder().alleles(Collections.singletonList(Allele.SV_SIMPLE_DEL)), true},
                {new GenotypeBuilder().alleles(Lists.newArrayList(Allele.REF_N, Allele.REF_N)), false},
                {new GenotypeBuilder().alleles(Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL)), true},
                {new GenotypeBuilder().alleles(Lists.newArrayList(Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DEL)), true},
                {new GenotypeBuilder().alleles(Lists.newArrayList(Allele.NO_CALL, Allele.SV_SIMPLE_DEL)), true},
                {new GenotypeBuilder().alleles(Lists.newArrayList(null, Allele.SV_SIMPLE_DEL)), true},
        };
    }

    @Test(dataProvider= "testIsAltGenotypeData")
    public void testIsAltGenotype(final GenotypeBuilder genotype, final boolean expected) {
        Assert.assertEquals(SVCallRecordUtils.containsAltAllele(genotype.make()), expected);
    }

    @DataProvider(name = "testSortAllelesData")
    public Object[][] testSortAllelesData() {
        return new Object[][]{
                {Collections.emptyList(), Collections.emptyList()},
                {Collections.singletonList(null), Collections.singletonList(null)},
                {Collections.singletonList(Allele.SV_SIMPLE_DEL), Collections.singletonList(Allele.SV_SIMPLE_DEL)},
                {Lists.newArrayList(null, Allele.SV_SIMPLE_DEL), Lists.newArrayList(null, Allele.SV_SIMPLE_DEL)},
                {Lists.newArrayList(Allele.SV_SIMPLE_DEL, null), Lists.newArrayList(null, Allele.SV_SIMPLE_DEL)},
                {Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL), Lists.newArrayList(Allele.SV_SIMPLE_DEL, Allele.REF_N)},
                {Lists.newArrayList(Allele.SV_SIMPLE_DEL, Allele.REF_N), Lists.newArrayList(Allele.SV_SIMPLE_DEL, Allele.REF_N)},
                {Lists.newArrayList(Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DUP), Lists.newArrayList(Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DUP)},
                {Lists.newArrayList(Allele.SV_SIMPLE_DUP, Allele.SV_SIMPLE_DEL), Lists.newArrayList(Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DUP)},
        };
    }

    @Test(dataProvider= "testSortAllelesData")
    public void testSortAlleles(final Collection<Allele> alleles, final List<Allele> expected) {
        Assert.assertEquals(SVCallRecordUtils.sortAlleles(alleles), expected);
    }



    @DataProvider(name = "testCreateData")
    public Object[][] testCreateData() {
        return new Object[][]{
                {
                    SVTestUtils.newVariantContext("var1", "chr1", 1000, 1999,
                            ALLELES_DEL, Lists.newArrayList(GENOTYPE_DEL_1, GENOTYPE_DEL_2), 1000, "+-",
                            StructuralVariantType.DEL, Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
                            null, null, Collections.singletonMap(GATKSVVCFConstants.COPY_NUMBER_FORMAT, 2)),
                        false,
                        new SVCallRecord("var1", "chr1", 1000, true, "chr1", 1999, false, StructuralVariantType.DEL, 1000,
                            Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM), ALLELES_DEL, Lists.newArrayList(GENOTYPE_DEL_1, GENOTYPE_DEL_2),
                                Collections.emptyMap())
                },
                {
                        SVTestUtils.newVariantContext("var1", "chr1", 1000, 1999,
                                ALLELES_DEL, Lists.newArrayList(GENOTYPE_DEL_1, GENOTYPE_DEL_2), 1000, "+-",
                                StructuralVariantType.DEL, Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
                                null, null, Collections.singletonMap(GATKSVVCFConstants.COPY_NUMBER_FORMAT, 2)),
                        true,
                        new SVCallRecord("var1", "chr1", 1000, true, "chr1", 1999, false, StructuralVariantType.DEL, 1000,
                                Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM), ALLELES_DEL, Lists.newArrayList(GENOTYPE_DEL_1, GENOTYPE_DEL_2),
                                Collections.singletonMap(GATKSVVCFConstants.COPY_NUMBER_FORMAT, 2))
                },
                {
                        SVTestUtils.newVariantContext("var2", "chr1", 1000, 1000,
                                ALLELES_INS, Lists.newArrayList(GENOTYPE_INS_1, GENOTYPE_INS_2), 500, "+-",
                                StructuralVariantType.INS, Collections.singletonList("pesr"),
                                null, null, Collections.emptyMap()),
                        false,
                        new SVCallRecord("var2", "chr1", 1000, true, "chr1", 1000, false, StructuralVariantType.INS, 500,
                                Collections.singletonList("pesr"), ALLELES_INS, Lists.newArrayList(GENOTYPE_INS_1, GENOTYPE_INS_2),
                                Collections.emptyMap())
                },
                {
                        SVTestUtils.newVariantContext("var3", "chr1", 1000, 1000,
                                ALLELES_BND, Collections.singletonList(GENOTYPE_BND_1), null, "++",
                                StructuralVariantType.BND, Collections.singletonList("pesr"),
                                "chrX", 2000, Collections.emptyMap()),
                        false,
                        new SVCallRecord("var3", "chr1", 1000, true, "chrX", 2000, true, StructuralVariantType.BND, null,
                                Collections.singletonList("pesr"), ALLELES_BND, Collections.singletonList(GENOTYPE_BND_1),
                                Collections.emptyMap())
                },
        };
    }

    @Test(dataProvider= "testCreateData")
    public void testCreate(final VariantContext variant, final boolean keepAttr, final SVCallRecord expected) {
        final List<String> excludedAttributes = new ArrayList<>(variant.getAttributes().keySet());
        excludedAttributes.remove(GATKSVVCFConstants.COPY_NUMBER_FORMAT); // Remove any attributes set in test data
        SVTestUtils.assertEqualsExceptExcludedAttributes(SVCallRecordUtils.create(variant, keepAttr), expected, excludedAttributes);
        if (keepAttr) {
            SVTestUtils.assertEqualsExceptExcludedAttributes(SVCallRecordUtils.create(variant), expected, excludedAttributes);
        }
    }

    @Test
    public void testCreateKeepAttr() {
        final VariantContext variant = SVTestUtils.newVariantContext("var1", "chr1", 1000, 1999,
                ALLELES_DEL, Lists.newArrayList(GENOTYPE_DEL_1, GENOTYPE_DEL_2), 1000, "+-",
                StructuralVariantType.DEL, Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
                null, null, Collections.singletonMap(GATKSVVCFConstants.COPY_NUMBER_FORMAT, 2));
        final SVCallRecord record = new SVCallRecord("var1", "chr1", 1000, true, "chr1", 1999, false, StructuralVariantType.DEL, 1000,
                Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM), ALLELES_DEL, Lists.newArrayList(GENOTYPE_DEL_1, GENOTYPE_DEL_2),
                Collections.singletonMap(GATKSVVCFConstants.COPY_NUMBER_FORMAT, 2));
        final List<String> excludedAttributes = new ArrayList<>(variant.getAttributes().keySet());
        excludedAttributes.remove(GATKSVVCFConstants.COPY_NUMBER_FORMAT);
        SVTestUtils.assertEqualsExceptExcludedAttributes(SVCallRecordUtils.create(variant, true), record, excludedAttributes);
        SVTestUtils.assertEqualsExceptExcludedAttributes(SVCallRecordUtils.create(variant), record, excludedAttributes);
    }
}