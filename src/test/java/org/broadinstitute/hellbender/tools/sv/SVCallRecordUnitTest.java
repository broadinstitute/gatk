package org.broadinstitute.hellbender.tools.sv;

import com.google.common.collect.Lists;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.GenotypesContext;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Collections;

public class SVCallRecordUnitTest {
    @DataProvider(name = "testIsIntrachomosomalData")
    public Object[][] testIsIntrachomosomalData() {
        return new Object[][]{
                {SVTestUtils.newCallRecordWithCoordinates("", "chr1", 1000, "chr1", 1999), true},
                {SVTestUtils.newCallRecordWithCoordinates("", "chr1", 1000, "chrX", 1999), false},
                {SVTestUtils.newCallRecordWithCoordinates("", "chrX", 1000, "chr1", 1999), false}
        };
    }

    @Test(dataProvider= "testIsIntrachomosomalData")
    public void testIsIntrachomosomal(final SVCallRecord record, final boolean expected) {
        Assert.assertEquals(record.isIntrachromosomal(), expected);
    }

    @DataProvider(name = "testIsDepthOnlyData")
    public Object[][] testIsDepthOnlyData() {
        return new Object[][]{
                {SVTestUtils.newCallRecordWithAlgorithms(Collections.emptyList()), false},
                {SVTestUtils.newCallRecordWithAlgorithms(Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM)), true},
                {SVTestUtils.newCallRecordWithAlgorithms(SVTestUtils.PESR_ONLY_ALGORITHM_LIST), false},
                {SVTestUtils.newCallRecordWithAlgorithms(Lists.newArrayList(GATKSVVCFConstants.DEPTH_ALGORITHM, SVTestUtils.PESR_ALGORITHM)), false}
        };
    }

    @Test(dataProvider= "testIsDepthOnlyData")
    public void testIsDepthOnly(final SVCallRecord record, final boolean expected) {
        Assert.assertEquals(record.isDepthOnly(), expected);
    }

    @DataProvider(name = "testIsCNVData")
    public Object[][] testIsCNVData() {
        return new Object[][]{
                {SVTestUtils.newPESRCallRecordWithIntervalAndType(1000, 1999, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL), true},
                {SVTestUtils.newPESRCallRecordWithIntervalAndType(1000, 1999, GATKSVVCFConstants.StructuralVariantAnnotationType.DUP), true},
                {SVTestUtils.newPESRCallRecordWithIntervalAndType(1000, 1999, GATKSVVCFConstants.StructuralVariantAnnotationType.CNV), true},
                {SVTestUtils.newPESRCallRecordWithIntervalAndType(1000, 1999, GATKSVVCFConstants.StructuralVariantAnnotationType.INV), false},
                {SVTestUtils.newPESRCallRecordWithIntervalAndType(1000, 1000, GATKSVVCFConstants.StructuralVariantAnnotationType.BND), false},
                {SVTestUtils.newPESRCallRecordWithIntervalAndType(1000, 1000, GATKSVVCFConstants.StructuralVariantAnnotationType.INS), false},
                {SVTestUtils.newPESRCallRecordWithIntervalAndType(1000, 1000, GATKSVVCFConstants.StructuralVariantAnnotationType.CPX), false},
                {SVTestUtils.newPESRCallRecordWithIntervalAndType(1000, 1000, GATKSVVCFConstants.StructuralVariantAnnotationType.CTX), false}
        };
    }

    @Test(dataProvider= "testIsCNVData")
    public void testIsCNV(final SVCallRecord record, final boolean expected) {
        Assert.assertEquals(record.isSimpleCNV(), expected);
    }

    @DataProvider(name = "testCreateInvalidCoordinatesData")
    public Object[][] testCreateInvalidCoordinatesData() {
        return new Object[][]{
                {"chr1", 0, "chr1", 248956422},
                {"chr1", 1, "chr1", 248956423},
                {"chr1", 1, "chr1", 248956423},
                {"chr1", 1, "chr2", 242193530},
                {"chr1", 2, "chr1", 1},
                {"chr2", 1, "chr1", 2}
        };
    }

    @Test(dataProvider="testCreateInvalidCoordinatesData", expectedExceptions = { IllegalArgumentException.class })
    public void testCreateInvalidCoordinates(final String contigA, final int posA, final String contigB, final int posB) {
        new SVCallRecord("var1", contigA, posA, true, contigB, posB, false, GATKSVVCFConstants.StructuralVariantAnnotationType.BND,
                null, Collections.emptyList(), null, SVTestUtils.PESR_ONLY_ALGORITHM_LIST, Collections.emptyList(), Collections.emptyList(),
                Collections.emptyMap(), Collections.emptySet(), null, SVTestUtils.hg38Dict);
        Assert.fail("Expected exception not thrown");
    }

    @DataProvider(name = "testCreateValidCoordinatesData")
    public Object[][] testCreateValidCoordinatesData() {
        return new Object[][]{
                {"chr1", 1, "chr1", 1},  // Start == END should be valid, e.g. for insertions
                {"chr1", 1, "chr1", 2},
                {"chr1", 2, "chr2", 1}
        };
    }

    @Test(dataProvider="testCreateValidCoordinatesData")
    public void testCreateValidCoordinates(final String contigA, final int posA, final String contigB, final int posB) {
        new SVCallRecord("var1", contigA, posA, true, contigB, posB, false, GATKSVVCFConstants.StructuralVariantAnnotationType.BND,
                null, Collections.emptyList(), null, SVTestUtils.PESR_ONLY_ALGORITHM_LIST, Collections.emptyList(), Collections.emptyList(),
                Collections.emptyMap(), Collections.emptySet(), null, SVTestUtils.hg38Dict);
    }

    @Test
    public void testGetters() {
        final SVCallRecord record = new SVCallRecord("var1", "chr1", 100, true, "chr1", 200, false, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL,
                GATKSVVCFConstants.ComplexVariantSubtype.dDUP, Lists.newArrayList(new SVCallRecord.ComplexEventInterval("DUP_chr1:100-200")), null, SVTestUtils.PESR_ONLY_ALGORITHM_LIST, Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL),
                GenotypesContext.create(GenotypeBuilder.create("sample1", Lists.newArrayList(Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DEL))),
                Collections.singletonMap("TEST_KEY", "TEST_VALUE"), Collections.singleton("TEST_FILTER"), Double.valueOf(30), SVTestUtils.hg38Dict);
        Assert.assertEquals(record.getId(), "var1");
        Assert.assertEquals(record.getContigA(), "chr1");
        Assert.assertEquals(record.getPositionA(), 100);
        Assert.assertEquals(record.getStrandA(), Boolean.TRUE);
        Assert.assertEquals(record.getContigB(), "chr1");
        Assert.assertEquals(record.getPositionB(), 200);
        Assert.assertEquals(record.getStrandB(), Boolean.FALSE);
        Assert.assertEquals(record.getAlgorithms(), SVTestUtils.PESR_ONLY_ALGORITHM_LIST);
        Assert.assertEquals(record.getGenotypes().get("sample1").getAlleles(), Lists.newArrayList(Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DEL));
        Assert.assertEquals(record.getAttributes(), Collections.singletonMap("TEST_KEY", "TEST_VALUE"));
        Assert.assertEquals(record.getAlleles(), Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL));
        Assert.assertEquals(record.getFilters(), Collections.singleton("TEST_FILTER"));
        Assert.assertEquals(record.getLog10PError(), Double.valueOf(30));
        Assert.assertEquals(record.getComplexSubtype(), GATKSVVCFConstants.ComplexVariantSubtype.dDUP);
        Assert.assertEquals(record.getComplexEventIntervals(), Lists.newArrayList(new SVCallRecord.ComplexEventInterval("DUP_chr1:100-200")));
    }

    @Test
    public void testComplexEventInterval() {
        final SVCallRecord.ComplexEventInterval cpx1 = new SVCallRecord.ComplexEventInterval("DEL_chr1:100-200");
        Assert.assertEquals(cpx1.getIntervalType(), GATKSVVCFConstants.StructuralVariantAnnotationType.DEL);
        Assert.assertEquals(cpx1.getInterval(), new SimpleInterval("chr1", 100, 200));
        Assert.assertEquals(cpx1.getContig(), "chr1");
        Assert.assertEquals(cpx1.getStart(), 100);
        Assert.assertEquals(cpx1.getEnd(), 200);
        Assert.assertEquals(cpx1.encode(), "DEL_chr1:100-200");
    }

    @DataProvider(name = "testComplexEventIntervalEqualsData")
    public Object[][] testComplexEventIntervalEqualsData() {
        return new Object[][]{
                {new SVCallRecord.ComplexEventInterval("DEL_chr1:100-200"), true},
                {new SVCallRecord.ComplexEventInterval(GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, new SimpleInterval("chr1", 100, 200)), true},
                {new SVCallRecord.ComplexEventInterval(GATKSVVCFConstants.StructuralVariantAnnotationType.DUP, new SimpleInterval("chr1", 100, 200)), false},
                {new SVCallRecord.ComplexEventInterval(GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, new SimpleInterval("chr2", 100, 200)), false},
                {new SVCallRecord.ComplexEventInterval(GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, new SimpleInterval("chr1", 101, 200)), false},
                {new SVCallRecord.ComplexEventInterval(GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, new SimpleInterval("chr1", 100, 201)), false}
        };
    }
    @Test(dataProvider="testComplexEventIntervalEqualsData")
    public void testComplexEventIntervalEquals(final SVCallRecord.ComplexEventInterval cpx2, final boolean expected) {
        final SVCallRecord.ComplexEventInterval cpx1 = new SVCallRecord.ComplexEventInterval("DEL_chr1:100-200");
        Assert.assertEquals(cpx1.equals(cpx2), expected);
        if (expected) {
            Assert.assertEquals(cpx1.hashCode(), cpx2.hashCode());
        }
    }
}