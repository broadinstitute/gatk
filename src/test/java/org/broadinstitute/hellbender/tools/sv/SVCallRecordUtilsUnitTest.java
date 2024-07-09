package org.broadinstitute.hellbender.tools.sv;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFHeaderLines;
import org.broadinstitute.hellbender.tools.sv.cluster.PloidyTable;
import org.broadinstitute.hellbender.utils.variant.GATKSVVariantContextUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;
import java.util.stream.Collectors;

public class SVCallRecordUtilsUnitTest {

    private static final List<Allele> ALLELES_DEL = Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL);
    private static final List<Allele> ALLELES_INS = Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_INS);
    private static final List<Allele> ALLELES_BND = Lists.newArrayList(Allele.REF_N, GATKSVVariantContextUtils.BND_ALLELE);
    private static final List<Allele> ALLELES_CTX = Lists.newArrayList(Allele.REF_N, GATKSVVariantContextUtils.CTX_ALLELE);
    private static final List<Allele> ALLELES_CPX = Lists.newArrayList(Allele.REF_N, GATKSVVariantContextUtils.CPX_ALLELE);

    private static final Map<String, Object> TEST_ATTRIBUTES = Collections.singletonMap("TEST_KEY", "TEST_VAL");
    private static final Map<String, Object> TEST_ATTRIBUTES_CPX = Lists.newArrayList(
            new AbstractMap.SimpleImmutableEntry<String, Object>("TEST_KEY", "TEST_VAL"),
            new AbstractMap.SimpleImmutableEntry<String, Object>(GATKSVVCFConstants.CPX_TYPE, GATKSVVCFConstants.ComplexVariantSubtype.dDUP.toString())
            ).stream().collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue));
    private static final Map<String, Object> TEST_ATTRIBUTES_CTX = Map.of(
            "TEST_KEY", "TEST_VAL",
            GATKSVVCFConstants.CPX_TYPE, "CTX_PP/QQ"
    );

    private static final Genotype GENOTYPE_DEL_1 = new GenotypeBuilder("sample1")
            .alleles(Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL))
            .attribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT, 1)
            .attribute(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, 2)
            .make();
    private static final Genotype GENOTYPE_DEL_2 = new GenotypeBuilder("sample2")
            .alleles(Lists.newArrayList(Allele.REF_N, Allele.REF_N))
            .attribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT, 2)
            .attribute(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, 2)
            .make();
    private static final Genotype GENOTYPE_DEL_3 = new GenotypeBuilder("sample3")
            .alleles(Lists.newArrayList(Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DEL))
            .attribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT, 0)
            .attribute(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, 2)
            .make();

    private static final Genotype GENOTYPE_INS_1 = new GenotypeBuilder("sample1")
            .alleles(Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_INS)).make();
    private static final Genotype GENOTYPE_INS_2 = new GenotypeBuilder("sample2")
            .alleles(Lists.newArrayList(Allele.SV_SIMPLE_INS, Allele.SV_SIMPLE_INS)).make();
    private static final Genotype GENOTYPE_BND_1 = new GenotypeBuilder("sample1")
            .alleles(Lists.newArrayList(Allele.REF_N, GATKSVVariantContextUtils.BND_ALLELE)).make();
    private static final Genotype GENOTYPE_CTX_1 = new GenotypeBuilder("sample1")
            .alleles(Lists.newArrayList(Allele.REF_N, GATKSVVariantContextUtils.CTX_ALLELE)).make();
    private static final Genotype GENOTYPE_CPX_1 = new GenotypeBuilder("sample1")
            .alleles(Lists.newArrayList(Allele.REF_N, GATKSVVariantContextUtils.CPX_ALLELE)).make();

    private static final Comparator<SVCallRecord> RECORD_COMPARATOR = SVCallRecordUtils.getCallComparator(SVTestUtils.hg38Dict);

    @DataProvider(name = "testGetVariantBuilderData")
    public Object[][] testGetVariantBuilderData() {
        return new Object[][]{
                // DEL
                {
                        new SVCallRecord("var1", "chr1", 1000, true, "chr1", 1999, false, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, null, Collections.emptyList(), 1000,
                                SVTestUtils.DEPTH_ONLY_ALGORITHM_LIST,
                                ALLELES_DEL,
                                Lists.newArrayList(GENOTYPE_DEL_1, GENOTYPE_DEL_2), Collections.emptyMap(), Collections.singleton("TEST_FILTER"), Double.valueOf(-3)),
                        new VariantContextBuilder("", "chr1", 1000, 1999, ALLELES_DEL)
                                .id("var1")
                                .genotypes(GENOTYPE_DEL_1, GENOTYPE_DEL_2)
                                .attribute(VCFConstants.END_KEY, 1999)
                                .attribute(GATKSVVCFConstants.ALGORITHMS_ATTRIBUTE, SVTestUtils.DEPTH_ONLY_ALGORITHM_LIST)
                                .attribute(GATKSVVCFConstants.SVTYPE, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL)
                                .filter("TEST_FILTER")
                                .log10PError(Double.valueOf(-3))
                                .make(),
                        Collections.emptyList()
                },
                // DEL w/ null ref allele
                {
                        new SVCallRecord("var1", "chr1", 1000, true, "chr1", 1999, false, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, null, Collections.emptyList(), 1000,
                                SVTestUtils.DEPTH_ONLY_ALGORITHM_LIST,
                                Collections.singletonList(Allele.SV_SIMPLE_DEL),
                                Collections.singletonList(GENOTYPE_DEL_3),
                                Collections.emptyMap(), Collections.emptySet(), null),
                        new VariantContextBuilder("", "chr1", 1000, 1999, ALLELES_DEL)
                                .id("var1")
                                .genotypes(GENOTYPE_DEL_3)
                                .attribute(VCFConstants.END_KEY, 1999)
                                .attribute(GATKSVVCFConstants.ALGORITHMS_ATTRIBUTE, SVTestUtils.DEPTH_ONLY_ALGORITHM_LIST)
                                .attribute(GATKSVVCFConstants.SVTYPE, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL)
                                .make(),
                        Collections.emptyList()
                },
                // INS
                {
                        new SVCallRecord("var2", "chr1", 1000, true, "chr1", 1000, false, GATKSVVCFConstants.StructuralVariantAnnotationType.INS, null, Collections.emptyList(), 500,
                                SVTestUtils.PESR_ONLY_ALGORITHM_LIST,
                                ALLELES_INS,
                                Lists.newArrayList(GENOTYPE_INS_1),
                                Collections.emptyMap(), Collections.emptySet(), null),
                        new VariantContextBuilder("", "chr1", 1000, 1000, ALLELES_INS)
                                .id("var2")
                                .genotypes(GENOTYPE_INS_1)
                                .attribute(VCFConstants.END_KEY, 1000)
                                .attribute(GATKSVVCFConstants.STRANDS_ATTRIBUTE, "+-")
                                .attribute(GATKSVVCFConstants.ALGORITHMS_ATTRIBUTE, SVTestUtils.PESR_ONLY_ALGORITHM_LIST)
                                .attribute(GATKSVVCFConstants.SVTYPE, GATKSVVCFConstants.StructuralVariantAnnotationType.INS)
                                .make(),
                        Collections.emptyList()
                },
                // INS, flipped strands
                {
                        new SVCallRecord("var2", "chr1", 1000, false, "chr1", 1000, true, GATKSVVCFConstants.StructuralVariantAnnotationType.INS, null, Collections.emptyList(), 500,
                                SVTestUtils.PESR_ONLY_ALGORITHM_LIST,
                                ALLELES_INS,
                                Lists.newArrayList(GENOTYPE_INS_1),
                                Collections.emptyMap(), Collections.emptySet(), null),
                        new VariantContextBuilder("", "chr1", 1000, 1000, ALLELES_INS)
                                .id("var2")
                                .genotypes(GENOTYPE_INS_1)
                                .attribute(VCFConstants.END_KEY, 1000)
                                .attribute(GATKSVVCFConstants.STRANDS_ATTRIBUTE, "-+")
                                .attribute(GATKSVVCFConstants.ALGORITHMS_ATTRIBUTE, SVTestUtils.PESR_ONLY_ALGORITHM_LIST)
                                .attribute(GATKSVVCFConstants.SVTYPE, GATKSVVCFConstants.StructuralVariantAnnotationType.INS)
                                .make(),
                        Collections.emptyList()
                },
                // INS, undefined length
                {
                        new SVCallRecord("var2", "chr1", 1000, true, "chr1", 1000, false, GATKSVVCFConstants.StructuralVariantAnnotationType.INS, null, Collections.emptyList(), null,
                                SVTestUtils.PESR_ONLY_ALGORITHM_LIST,
                                ALLELES_INS,
                                Lists.newArrayList(GENOTYPE_INS_1),
                                Collections.emptyMap(), Collections.emptySet(), null),
                        new VariantContextBuilder("", "chr1", 1000, 1000, ALLELES_INS)
                                .id("var2")
                                .genotypes(GENOTYPE_INS_1)
                                .attribute(VCFConstants.END_KEY, 1000)
                                .attribute(GATKSVVCFConstants.STRANDS_ATTRIBUTE, "+-")
                                .attribute(GATKSVVCFConstants.ALGORITHMS_ATTRIBUTE, SVTestUtils.PESR_ONLY_ALGORITHM_LIST)
                                .attribute(GATKSVVCFConstants.SVTYPE, GATKSVVCFConstants.StructuralVariantAnnotationType.INS)
                                .make(),
                        Collections.emptyList()
                },
                // BND
                {
                        new SVCallRecord("var_bnd", "chr1", 1000, false, "chr2", 1999, true, GATKSVVCFConstants.StructuralVariantAnnotationType.BND, null, Collections.emptyList(), null,
                                SVTestUtils.PESR_ONLY_ALGORITHM_LIST,
                                ALLELES_BND,
                                Lists.newArrayList(GENOTYPE_BND_1),
                                Collections.emptyMap(), Collections.emptySet(), null),
                        new VariantContextBuilder("", "chr1", 1000, 1000, ALLELES_BND)
                                .id("var_bnd")
                                .genotypes(GENOTYPE_BND_1)
                                .attribute(VCFConstants.END_KEY, 1000)
                                .attribute(GATKSVVCFConstants.ALGORITHMS_ATTRIBUTE, SVTestUtils.PESR_ONLY_ALGORITHM_LIST)
                                .attribute(GATKSVVCFConstants.STRANDS_ATTRIBUTE, "-+")
                                .attribute(GATKSVVCFConstants.SVTYPE, GATKSVVCFConstants.StructuralVariantAnnotationType.BND)
                                .attribute(GATKSVVCFConstants.CONTIG2_ATTRIBUTE, "chr2")
                                .attribute(GATKSVVCFConstants.END2_ATTRIBUTE, 1999)
                                .make(),
                        Collections.emptyList()
                },
                // CTX
                {
                        new SVCallRecord("var_ctx", "chr1", 1000, false, "chr2", 1999, true, GATKSVVCFConstants.StructuralVariantAnnotationType.CTX, null, Collections.emptyList(), null,
                                SVTestUtils.PESR_ONLY_ALGORITHM_LIST,
                                ALLELES_CTX,
                                Lists.newArrayList(GENOTYPE_CTX_1),
                                Collections.emptyMap(), Collections.emptySet(), null),
                        new VariantContextBuilder("", "chr1", 1000, 1000, ALLELES_CTX)
                                .id("var_ctx")
                                .genotypes(GENOTYPE_CTX_1)
                                .attribute(VCFConstants.END_KEY, 1000)
                                .attribute(GATKSVVCFConstants.ALGORITHMS_ATTRIBUTE, SVTestUtils.PESR_ONLY_ALGORITHM_LIST)
                                .attribute(GATKSVVCFConstants.STRANDS_ATTRIBUTE, "-+")
                                .attribute(GATKSVVCFConstants.SVTYPE, GATKSVVCFConstants.StructuralVariantAnnotationType.CTX)
                                .attribute(GATKSVVCFConstants.CONTIG2_ATTRIBUTE, "chr2")
                                .attribute(GATKSVVCFConstants.END2_ATTRIBUTE, 1999)
                                .make(),
                        Collections.emptyList()
                },
                // CPX
                {
                        new SVCallRecord("var_cpx", "chr1", 1000, null, "chr1", 1000, null,
                                GATKSVVCFConstants.StructuralVariantAnnotationType.CPX,
                                GATKSVVCFConstants.ComplexVariantSubtype.dDUP_iDEL,
                                Lists.newArrayList(SVCallRecord.ComplexEventInterval.decode("DUP_chr1:5000-5100", SVTestUtils.hg38Dict), SVCallRecord.ComplexEventInterval.decode("DEL_chr2:100-200", SVTestUtils.hg38Dict)),
                                100,
                                SVTestUtils.PESR_ONLY_ALGORITHM_LIST,
                                ALLELES_CPX,
                                Lists.newArrayList(GENOTYPE_CPX_1),
                                Collections.emptyMap(), Collections.emptySet(), null),
                        new VariantContextBuilder("", "chr1", 1000, 1000, ALLELES_CPX)
                                .id("var_cpx")
                                .genotypes(GENOTYPE_CPX_1)
                                .attribute(VCFConstants.END_KEY, 1000)
                                .attribute(GATKSVVCFConstants.ALGORITHMS_ATTRIBUTE, SVTestUtils.PESR_ONLY_ALGORITHM_LIST)
                                .attribute(GATKSVVCFConstants.SVTYPE, GATKSVVCFConstants.StructuralVariantAnnotationType.CPX)
                                .attribute(GATKSVVCFConstants.CPX_TYPE, GATKSVVCFConstants.ComplexVariantSubtype.dDUP_iDEL.name())
                                .attribute(GATKSVVCFConstants.SVLEN, 100)
                                .attribute(GATKSVVCFConstants.CPX_INTERVALS, "DEL_chr2:100-200,DUP_chr1:5000-5100")
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
    public void testGetVariantBuilderHasSanitizedNullAttributes() {
        final SVCallRecord record = new SVCallRecord("var3", "chr1", 1000, false, "chr2", 1999, true, GATKSVVCFConstants.StructuralVariantAnnotationType.BND, null, Collections.emptyList(), null,
                SVTestUtils.PESR_ONLY_ALGORITHM_LIST,
                ALLELES_BND,
                Lists.newArrayList(GENOTYPE_BND_1),
                Collections.emptyMap(), Collections.emptySet(), null);
        final VariantContext result = SVCallRecordUtils.getVariantBuilder(record).make();
        // BNDs shouldn't have a length
        Assert.assertFalse(result.hasAttribute(GATKSVVCFConstants.SVLEN));
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

        final Map<String, Integer> sample1PloidyMap = Collections.singletonMap("chr20", 2);
        final Map<String, Integer> sample2PloidyMap = Collections.singletonMap("chr20", 2);
        final Map<String, Integer> sample3PloidyMap = Collections.singletonMap("chr20", 2);
        final Map<String, Integer> sample4PloidyMap = Collections.singletonMap("chr20", 3);
        final Map<String, Map<String, Integer>> rawPloidyMap = new HashMap<>();
        rawPloidyMap.put("sample1", sample1PloidyMap);
        rawPloidyMap.put("sample2", sample2PloidyMap);
        rawPloidyMap.put("sample3", sample3PloidyMap);
        rawPloidyMap.put("sample4", sample4PloidyMap);
        final PloidyTable ploidyTable = new PloidyTable(rawPloidyMap);

        final VCFHeader headerWithCopyNumber = new VCFHeader();
        headerWithCopyNumber.addMetaDataLine(GATKSVVCFHeaderLines.getFormatLine(GATKSVVCFConstants.COPY_NUMBER_FORMAT));

        final SVCallRecord record1 = SVTestUtils.makeRecord("record1", "chr20", 1000, true,
                "chr20", 2000, false, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, null,
                Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM), Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL),
                Lists.newArrayList(new GenotypeBuilder(GENOTYPE_DEL_1), new GenotypeBuilder(GENOTYPE_DEL_2)));
        final GenotypesContext resultNoMissing = SVCallRecordUtils.populateGenotypesForMissingSamplesWithAlleles(
                record1, samplesNoMissing, true, ploidyTable, headerWithCopyNumber);
        Assert.assertEquals(resultNoMissing.size(), 2);
        for (final Genotype g : genotypesContext) {
            VariantContextTestUtils.assertGenotypesAreEqual(g, genotypesContext.get(g.getSampleName()));
        }

        final GenotypesContext result = SVCallRecordUtils.populateGenotypesForMissingSamplesWithAlleles(record1, samples, true, ploidyTable, headerWithCopyNumber);
        Assert.assertEquals(result.size(), 4);

        final Genotype g1 = result.get(GENOTYPE_DEL_1.getSampleName());
        final Genotype g2 = result.get(GENOTYPE_DEL_2.getSampleName());
        final Genotype g3 = result.get("sample3");
        final Genotype g4 = result.get("sample4");

        final List<Allele> alleles3 = Collections.nCopies(2, record1.getRefAllele());
        final Genotype g3Expected = new GenotypeBuilder("sample3", alleles3)
                .attribute(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, 2)
                .attribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT, 2)
                .make();
        final List<Allele> alleles4 = Collections.nCopies(3, record1.getRefAllele());
        final Genotype g4Expected = new GenotypeBuilder("sample4", alleles4)
                .attribute(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, 3)
                .attribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT, 3)
                .make();

        VariantContextTestUtils.assertGenotypesAreEqual(g1, GENOTYPE_DEL_1);
        VariantContextTestUtils.assertGenotypesAreEqual(g2, GENOTYPE_DEL_2);
        VariantContextTestUtils.assertGenotypesAreEqual(g3, g3Expected);
        VariantContextTestUtils.assertGenotypesAreEqual(g4, g4Expected);

        // Should omit CN from new genotypes since it's not in the given header
        final GenotypesContext resultNoCopyNumber = SVCallRecordUtils.populateGenotypesForMissingSamplesWithAlleles(record1, samples, true, ploidyTable, new VCFHeader());
        Assert.assertFalse(resultNoCopyNumber.get("sample3").hasExtendedAttribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT));
        Assert.assertFalse(resultNoCopyNumber.get("sample4").hasExtendedAttribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT));
    }

    @Test
    public void testCopyCallWithNewGenotypes() {
        final SVCallRecord record = new SVCallRecord("var1", "chr1", 1000, true, "chr1", 1999, false, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, null, Collections.emptyList(), 1000,
                SVTestUtils.DEPTH_ONLY_ALGORITHM_LIST,
                ALLELES_DEL,
                Lists.newArrayList(GENOTYPE_DEL_1, GENOTYPE_DEL_2),
                Collections.singletonMap(GATKSVVCFConstants.CLUSTER_MEMBER_IDS_KEY, Collections.singletonList("sample")), Collections.emptySet(), null);
        final GenotypesContext genotypes = GenotypesContext.copy(Collections.singletonList(GENOTYPE_DEL_3));
        final SVCallRecord result = SVCallRecordUtils.copyCallWithNewGenotypes(record, genotypes);
        final SVCallRecord expected = new SVCallRecord("var1", "chr1", 1000, true, "chr1", 1999, false, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, null, Collections.emptyList(), 1000,
                SVTestUtils.DEPTH_ONLY_ALGORITHM_LIST,
                ALLELES_DEL,
                genotypes,
                Collections.singletonMap(GATKSVVCFConstants.CLUSTER_MEMBER_IDS_KEY, Collections.singletonList("sample")), Collections.emptySet(), null);
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
                        SVTestUtils.newCallRecordWithCoordinatesAndType("var2", "chrX", 1000, "chr1", 1999, GATKSVVCFConstants.StructuralVariantAnnotationType.BND),
                        -1
                },
                {
                        SVTestUtils.newCallRecordWithCoordinatesAndType("var1", "chrX", 1000, "chr1", 1999, GATKSVVCFConstants.StructuralVariantAnnotationType.BND),
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
                        SVTestUtils.newCallRecordWithCoordinatesAndType("var2", "chr1", 1000, "chrX", 1999, GATKSVVCFConstants.StructuralVariantAnnotationType.BND),
                        -1
                },
                {
                        SVTestUtils.newCallRecordWithCoordinatesAndType("var1", "chr1", 1000, "chrX", 1999, GATKSVVCFConstants.StructuralVariantAnnotationType.BND),
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
                        SVTestUtils.newPESRCallRecordWithIntervalAndType(1000, 1999, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL),
                        SVTestUtils.newPESRCallRecordWithIntervalAndType(1000, 1999, GATKSVVCFConstants.StructuralVariantAnnotationType.DUP),
                        GATKSVVCFConstants.StructuralVariantAnnotationType.DEL.compareTo(GATKSVVCFConstants.StructuralVariantAnnotationType.DUP)
                },
                {
                        SVTestUtils.newPESRCallRecordWithIntervalAndType(1000, 1999, GATKSVVCFConstants.StructuralVariantAnnotationType.DUP),
                        SVTestUtils.newPESRCallRecordWithIntervalAndType(1000, 1999, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL),
                        GATKSVVCFConstants.StructuralVariantAnnotationType.DUP.compareTo(GATKSVVCFConstants.StructuralVariantAnnotationType.DEL)
                },
                {
                        SVTestUtils.newPESRCallRecordWithIntervalAndType(1000, 1000, GATKSVVCFConstants.StructuralVariantAnnotationType.BND),
                        SVTestUtils.newPESRCallRecordWithIntervalAndType(1000, 1000, GATKSVVCFConstants.StructuralVariantAnnotationType.INS),
                        GATKSVVCFConstants.StructuralVariantAnnotationType.BND.compareTo(GATKSVVCFConstants.StructuralVariantAnnotationType.INS)
                }
        };
    }

    @Test(dataProvider= "testGetCallComparatorData")
    public void testGetCallComparator(final SVCallRecord record1, final SVCallRecord record2, final int expected) {
        Assert.assertEquals(RECORD_COMPARATOR.compare(record1, record2), expected);
    }

    @Test
    public void testConvertInversionsToBreakends() {
        final SVCallRecord nonInversion = SVTestUtils.newPESRCallRecordWithIntervalAndType(1000, 1999, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL);
        final List<SVCallRecord> nonInversionResult = SVCallRecordUtils.convertInversionsToBreakends(nonInversion, SVTestUtils.hg38Dict).collect(Collectors.toList());
        Assert.assertEquals(nonInversionResult.size(), 1);
        Assert.assertNotNull(nonInversionResult.get(0));
        SVTestUtils.assertEqualsExceptMembership(nonInversionResult.get(0), nonInversion);

        final SVCallRecord inversion = new SVCallRecord("", "chr1", 1000, true, "chr1", 1999, true, GATKSVVCFConstants.StructuralVariantAnnotationType.INV, null, Collections.emptyList(), 1000,
                SVTestUtils.PESR_ONLY_ALGORITHM_LIST,
                Collections.emptyList(),
                Collections.emptyList(),
                Collections.emptyMap(), Collections.emptySet(), null);
        final List<SVCallRecord> inversionResult = SVCallRecordUtils.convertInversionsToBreakends(inversion, SVTestUtils.hg38Dict).collect(Collectors.toList());
        Assert.assertEquals(inversionResult.size(), 2);

        final SVCallRecord bnd1 = inversionResult.get(0);
        final SVCallRecord bnd2 = inversionResult.get(1);
        Assert.assertNotNull(bnd1);
        Assert.assertNotNull(bnd2);
        Assert.assertEquals(bnd1.getType(), GATKSVVCFConstants.StructuralVariantAnnotationType.BND);
        Assert.assertEquals(bnd2.getType(), GATKSVVCFConstants.StructuralVariantAnnotationType.BND);
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
                            GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
                            null, null,
                            TEST_ATTRIBUTES, -90.),
                        new SVCallRecord("var1", "chr1", 1000, true, "chr1", 1999, false, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, null, Collections.emptyList(), 1000,
                            Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM), ALLELES_DEL, Lists.newArrayList(GENOTYPE_DEL_1, GENOTYPE_DEL_2),
                                TEST_ATTRIBUTES, Collections.emptySet(), -90.)
                },
                {
                        SVTestUtils.newVariantContext("var2", "chr1", 1000, 1999,
                                ALLELES_DEL, Lists.newArrayList(GENOTYPE_DEL_1, GENOTYPE_DEL_2), 1000, "+-",
                                GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
                                null, null,
                                TEST_ATTRIBUTES, null),
                        new SVCallRecord("var2", "chr1", 1000, true, "chr1", 1999, false, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, null, Collections.emptyList(), 1000,
                                Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM), ALLELES_DEL, Lists.newArrayList(GENOTYPE_DEL_1, GENOTYPE_DEL_2),
                                TEST_ATTRIBUTES, Collections.emptySet(), null)
                },
                {
                        SVTestUtils.newVariantContext("var3", "chr1", 1000, 1000,
                                ALLELES_INS, Lists.newArrayList(GENOTYPE_INS_1, GENOTYPE_INS_2), 500, "+-",
                                GATKSVVCFConstants.StructuralVariantAnnotationType.INS, SVTestUtils.PESR_ONLY_ALGORITHM_LIST,
                                null, null, TEST_ATTRIBUTES, null),
                        new SVCallRecord("var3", "chr1", 1000, true, "chr1", 1000, false, GATKSVVCFConstants.StructuralVariantAnnotationType.INS, null, Collections.emptyList(), 500,
                                SVTestUtils.PESR_ONLY_ALGORITHM_LIST, ALLELES_INS, Lists.newArrayList(GENOTYPE_INS_1, GENOTYPE_INS_2),
                                TEST_ATTRIBUTES, Collections.emptySet(), null)
                },
                {
                        SVTestUtils.newVariantContext("var4", "chr1", 1000, 1000,
                                ALLELES_INS, Lists.newArrayList(GENOTYPE_INS_1, GENOTYPE_INS_2), 500, "-+",
                                GATKSVVCFConstants.StructuralVariantAnnotationType.INS, SVTestUtils.PESR_ONLY_ALGORITHM_LIST,
                                null, null, TEST_ATTRIBUTES, null),
                        new SVCallRecord("var4", "chr1", 1000, false, "chr1", 1000, true, GATKSVVCFConstants.StructuralVariantAnnotationType.INS, null, Collections.emptyList(), 500,
                                SVTestUtils.PESR_ONLY_ALGORITHM_LIST, ALLELES_INS, Lists.newArrayList(GENOTYPE_INS_1, GENOTYPE_INS_2),
                                TEST_ATTRIBUTES, Collections.emptySet(), null)
                },
                {
                        SVTestUtils.newVariantContext("var4b", "chr1", 1000, 1000,
                                ALLELES_INS, Lists.newArrayList(GENOTYPE_INS_1, GENOTYPE_INS_2), -1, "-+",
                                GATKSVVCFConstants.StructuralVariantAnnotationType.INS, SVTestUtils.PESR_ONLY_ALGORITHM_LIST,
                                null, null, TEST_ATTRIBUTES, null),
                        new SVCallRecord("var4b", "chr1", 1000, false, "chr1", 1000, true, GATKSVVCFConstants.StructuralVariantAnnotationType.INS, null, Collections.emptyList(), null,
                                SVTestUtils.PESR_ONLY_ALGORITHM_LIST, ALLELES_INS, Lists.newArrayList(GENOTYPE_INS_1, GENOTYPE_INS_2),
                                TEST_ATTRIBUTES, Collections.emptySet(), null)
                },
                {
                        SVTestUtils.newVariantContext("var5", "chr1", 1000, 1000,
                                ALLELES_BND, Collections.singletonList(GENOTYPE_BND_1), null, "++",
                                GATKSVVCFConstants.StructuralVariantAnnotationType.BND, SVTestUtils.PESR_ONLY_ALGORITHM_LIST,
                                "chrX", 2000, TEST_ATTRIBUTES, null),
                        new SVCallRecord("var5", "chr1", 1000, true, "chrX", 2000, true, GATKSVVCFConstants.StructuralVariantAnnotationType.BND, null, Collections.emptyList(), null,
                                SVTestUtils.PESR_ONLY_ALGORITHM_LIST, ALLELES_BND, Collections.singletonList(GENOTYPE_BND_1),
                                TEST_ATTRIBUTES, Collections.emptySet(), null)
                },
                {
                        SVTestUtils.newVariantContext("var6", "chr1", 1000, 1000,
                                ALLELES_BND, Collections.singletonList(GENOTYPE_BND_1), null, "++",
                                GATKSVVCFConstants.StructuralVariantAnnotationType.BND, SVTestUtils.PESR_ONLY_ALGORITHM_LIST,
                                "chrX", 2000, TEST_ATTRIBUTES, null),
                        new SVCallRecord("var6", "chr1", 1000, true, "chrX", 2000, true, GATKSVVCFConstants.StructuralVariantAnnotationType.BND, null, Collections.emptyList(), null,
                                SVTestUtils.PESR_ONLY_ALGORITHM_LIST, ALLELES_BND, Collections.singletonList(GENOTYPE_BND_1),
                                TEST_ATTRIBUTES, Collections.emptySet(), null)
                },
                {
                        SVTestUtils.newVariantContext("var7", "chr1", 1000, 1000,
                                ALLELES_CPX, Collections.singletonList(GENOTYPE_CPX_1), 250, null,
                                GATKSVVCFConstants.StructuralVariantAnnotationType.CPX, SVTestUtils.PESR_ONLY_ALGORITHM_LIST,
                                "chrX", 2000, TEST_ATTRIBUTES_CPX, null),
                        new SVCallRecord("var7", "chr1", 1000, null, "chr1", 1000, null, GATKSVVCFConstants.StructuralVariantAnnotationType.CPX, GATKSVVCFConstants.ComplexVariantSubtype.dDUP, Collections.emptyList(), 250,
                                SVTestUtils.PESR_ONLY_ALGORITHM_LIST, ALLELES_CPX, Collections.singletonList(GENOTYPE_CPX_1),
                                TEST_ATTRIBUTES, Collections.emptySet(), null)
                },
                {
                        SVTestUtils.newVariantContext("var8", "chr1", 1000, 2000,
                                ALLELES_CPX, Collections.singletonList(GENOTYPE_CPX_1), 250, null,
                                GATKSVVCFConstants.StructuralVariantAnnotationType.CPX, SVTestUtils.PESR_ONLY_ALGORITHM_LIST,
                                "chr1", null, TEST_ATTRIBUTES_CPX, null),
                        new SVCallRecord("var8", "chr1", 1000, null, "chr1", 2000, null, GATKSVVCFConstants.StructuralVariantAnnotationType.CPX, GATKSVVCFConstants.ComplexVariantSubtype.dDUP, Collections.emptyList(), 250,
                                SVTestUtils.PESR_ONLY_ALGORITHM_LIST, ALLELES_CPX, Collections.singletonList(GENOTYPE_CPX_1),
                                TEST_ATTRIBUTES, Collections.emptySet(), null)
                },
                {
                        SVTestUtils.newVariantContext("var9", "chr1", 1000, 2000,
                                ALLELES_CPX, Collections.singletonList(GENOTYPE_CPX_1), 250, null,
                                GATKSVVCFConstants.StructuralVariantAnnotationType.CPX, SVTestUtils.PESR_ONLY_ALGORITHM_LIST,
                                null, null, TEST_ATTRIBUTES_CPX, null),
                        new SVCallRecord("var9", "chr1", 1000, null, "chr1", 2000, null, GATKSVVCFConstants.StructuralVariantAnnotationType.CPX, GATKSVVCFConstants.ComplexVariantSubtype.dDUP, Collections.emptyList(), 250,
                                SVTestUtils.PESR_ONLY_ALGORITHM_LIST, ALLELES_CPX, Collections.singletonList(GENOTYPE_CPX_1),
                                TEST_ATTRIBUTES, Collections.emptySet(), null)
                },
                {
                        SVTestUtils.newVariantContext("var10", "chr1", 1000, 2000,
                                ALLELES_CPX, Collections.singletonList(GENOTYPE_CPX_1), 250, null,
                                GATKSVVCFConstants.StructuralVariantAnnotationType.CPX, SVTestUtils.PESR_ONLY_ALGORITHM_LIST,
                                "chrX", 2000,
                                Map.of("TEST_KEY", "TEST_VAL", GATKSVVCFConstants.CPX_TYPE, GATKSVVCFConstants.ComplexVariantSubtype.dDUP_iDEL.name(), GATKSVVCFConstants.CPX_INTERVALS, Arrays.asList("DUP_chr1:100-200", "DEL_chr2:300-400")),
                                null),
                        new SVCallRecord("var10", "chr1", 1000, null, "chr1", 2000, null,
                                GATKSVVCFConstants.StructuralVariantAnnotationType.CPX,
                                GATKSVVCFConstants.ComplexVariantSubtype.dDUP_iDEL,
                                Lists.newArrayList(SVCallRecord.ComplexEventInterval.decode("DUP_chr1:100-200", SVTestUtils.hg38Dict), SVCallRecord.ComplexEventInterval.decode("DEL_chr2:300-400", SVTestUtils.hg38Dict)),
                                250,
                                SVTestUtils.PESR_ONLY_ALGORITHM_LIST, ALLELES_CPX, Collections.singletonList(GENOTYPE_CPX_1),
                                TEST_ATTRIBUTES, Collections.emptySet(), null)
                },
                {
                        SVTestUtils.newVariantContext("var11", "chr1", 1000, 1000,
                                ALLELES_CTX, Collections.singletonList(GENOTYPE_CTX_1), null, "++",
                                GATKSVVCFConstants.StructuralVariantAnnotationType.CTX, SVTestUtils.PESR_ONLY_ALGORITHM_LIST,
                                "chrX", 2000, TEST_ATTRIBUTES_CTX, null),
                        new SVCallRecord("var11", "chr1", 1000, true, "chrX", 2000, true, GATKSVVCFConstants.StructuralVariantAnnotationType.CTX, GATKSVVCFConstants.ComplexVariantSubtype.CTX_PP_QQ, Collections.emptyList(), null,
                                SVTestUtils.PESR_ONLY_ALGORITHM_LIST, ALLELES_CTX, Collections.singletonList(GENOTYPE_CTX_1),
                                TEST_ATTRIBUTES, Collections.emptySet(), null)
                },
                // CPX with a null value in the CPX_INTERVALS list
                {
                        SVTestUtils.newVariantContext("var12", "chr1", 1000, 1000,
                                ALLELES_CTX, Collections.singletonList(GENOTYPE_CTX_1), null, "++",
                                GATKSVVCFConstants.StructuralVariantAnnotationType.CTX, SVTestUtils.PESR_ONLY_ALGORITHM_LIST,
                                "chrX", 2000,
                                Map.of("TEST_KEY", "TEST_VAL", GATKSVVCFConstants.CPX_TYPE, "CTX_PP/QQ", GATKSVVCFConstants.CPX_INTERVALS, Collections.emptyList()),
                                null),
                        new SVCallRecord("var12", "chr1", 1000, true, "chrX", 2000, true, GATKSVVCFConstants.StructuralVariantAnnotationType.CTX, GATKSVVCFConstants.ComplexVariantSubtype.CTX_PP_QQ, Collections.emptyList(), null,
                                SVTestUtils.PESR_ONLY_ALGORITHM_LIST, ALLELES_CTX, Collections.singletonList(GENOTYPE_CTX_1),
                                TEST_ATTRIBUTES, Collections.emptySet(), null)
                },
        };
    }

    @Test(dataProvider= "testCreateData")
    public void testCreate(final VariantContext variant, final SVCallRecord expected) {
        final List<String> excludedAttributes = new ArrayList<>(variant.getAttributes().keySet());
        final SVCallRecord resultDropAttr = SVCallRecordUtils.create(variant, false, SVTestUtils.hg38Dict);
        SVTestUtils.assertEqualsExceptExcludedAttributes(resultDropAttr, expected, excludedAttributes);
        Assert.assertTrue(resultDropAttr.getAttributes().isEmpty());

        final SVCallRecord resultKeepAttr = SVCallRecordUtils.create(variant, true, SVTestUtils.hg38Dict);
        SVTestUtils.assertEqualsExceptExcludedAttributes(resultKeepAttr, expected, Collections.emptyList());
    }

    @DataProvider(name = "testGetComplexSubtypeData")
    public Object[][] testGetComplexSubtypeData() {
        return new Object[][]{
                {new VariantContextBuilder()
                        .source("source")
                        .id("id")
                        .chr("chr1")
                        .start(2000)
                        .stop(3000)
                        .alleles(Arrays.asList(Allele.REF_N, Allele.create("<CPX>", false)))
                        .attributes(Map.of(
                                GATKSVVCFConstants.SVTYPE, GATKSVVCFConstants.StructuralVariantAnnotationType.CPX,
                                GATKSVVCFConstants.CPX_TYPE, "dupINVdup"
                        ))
                        .make(),
                        GATKSVVCFConstants.ComplexVariantSubtype.dupINVdup
                },
                {new VariantContextBuilder()
                        .source("source")
                        .id("id")
                        .chr("chr1")
                        .start(2000)
                        .stop(3000)
                        .alleles(Arrays.asList(Allele.REF_N, Allele.create("<CPX>", false)))
                        .attributes(Map.of(
                                GATKSVVCFConstants.SVTYPE, GATKSVVCFConstants.StructuralVariantAnnotationType.CPX,
                                GATKSVVCFConstants.CPX_TYPE, "CTX_PP/QQ"
                        ))
                        .make(),
                        GATKSVVCFConstants.ComplexVariantSubtype.CTX_PP_QQ
                },
                {new VariantContextBuilder()
                        .source("source")
                        .id("id")
                        .chr("chr1")
                        .start(2000)
                        .stop(3000)
                        .alleles(Arrays.asList(Allele.REF_N, Allele.create("<DEL>", false)))
                        .attributes(Map.of(
                                GATKSVVCFConstants.SVTYPE, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL
                        ))
                        .make(),
                        null
                }
        };
    }

    @Test(dataProvider= "testGetComplexSubtypeData")
    public void testGetComplexSubtype(final VariantContext variant, final GATKSVVCFConstants.ComplexVariantSubtype expected) {
        final GATKSVVCFConstants.ComplexVariantSubtype actual = SVCallRecordUtils.getComplexSubtype(variant);
        Assert.assertEquals(actual, expected);
    }

    @DataProvider(name = "testGetComplexSubtypeStringData")
    public Object[][] testGetComplexSubtypeStringData() {
        return new Object[][]{
                {GATKSVVCFConstants.ComplexVariantSubtype.CTX_PQ_QP, "CTX_PQ/QP"},
                {GATKSVVCFConstants.ComplexVariantSubtype.CTX_PP_QQ, "CTX_PP/QQ"},
                {GATKSVVCFConstants.ComplexVariantSubtype.INS_iDEL, "INS_iDEL"}
        };
    }

    @Test(dataProvider= "testGetComplexSubtypeStringData")
    public void testGetComplexSubtypeString(final GATKSVVCFConstants.ComplexVariantSubtype subtype, final String expected) {
        final String actual = SVCallRecordUtils.getComplexSubtypeString(subtype);
        Assert.assertEquals(actual, expected);
    }
}