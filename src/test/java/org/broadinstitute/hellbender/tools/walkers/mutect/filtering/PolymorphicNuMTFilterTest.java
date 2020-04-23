package org.broadinstitute.hellbender.tools.walkers.mutect.filtering;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;

import static org.broadinstitute.hellbender.GATKBaseTest.toolsTestDir;

public class PolymorphicNuMTFilterTest extends GATKBaseTest {
    private static File statsFile = new File(toolsTestDir, "mutect/mito/unfiltered.vcf.stats");
    private static VCFHeader header = VariantContextTestUtils.getVCFHeader(toolsTestDir + "mutect/mito/unfiltered.vcf");
    private static Mutect2FilteringEngine engine = new Mutect2FilteringEngine(new M2FiltersArgumentCollection(), header, statsFile);

    @Test(dataProvider = "getCoveragesAndExpected")
    public void testIsArtifact(final double medianAutosomalCoverage, final int altAlleleCount, final boolean expectedResult) {
        final PolymorphicNuMTFilter filter = new PolymorphicNuMTFilter(medianAutosomalCoverage);
        final VariantContext vc = VariantContextTestUtils.makeSomaticNuMTVariant(altAlleleCount);
        Assert.assertEquals(filter.isArtifact(vc, engine), expectedResult);
    }

    @DataProvider(name = "getCoveragesAndExpected")
    public Object[][] getCoveragesAndExpected() {
        return new Object[][] {
                {25, 11, true},
                {22, 11, true},
                {3, 11, false},
                {22, 46, true},
                {22, 47, false},
                {25, 51, true},
                {25, 52, false},
                {27, 55, true},
                {27, 56, false}
        };
    }
}