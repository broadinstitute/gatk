package org.broadinstitute.hellbender.tools.sv;

import com.google.common.collect.Lists;
import htsjdk.variant.variantcontext.StructuralVariantType;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
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
                {SVTestUtils.newCallRecordWithIntervalAndType(1000, 1999, StructuralVariantType.DEL), true},
                {SVTestUtils.newCallRecordWithIntervalAndType(1000, 1999, StructuralVariantType.DUP), true},
                {SVTestUtils.newCallRecordWithIntervalAndType(1000, 1999, StructuralVariantType.CNV), true},
                {SVTestUtils.newCallRecordWithIntervalAndType(1000, 1999, StructuralVariantType.INV), false},
                {SVTestUtils.newCallRecordWithIntervalAndType(1000, 1999, StructuralVariantType.BND), false},
                {SVTestUtils.newCallRecordWithIntervalAndType(1000, 1000, StructuralVariantType.INS), false}
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
        new SVCallRecord("var1", contigA, posA, true, contigB, posB, false, StructuralVariantType.BND,
                null, SVTestUtils.PESR_ONLY_ALGORITHM_LIST, Collections.emptyList(), Collections.emptyList(),
                Collections.emptyMap(), SVTestUtils.hg38Dict);
        Assert.fail("Expected exception not thrown");
    }
}