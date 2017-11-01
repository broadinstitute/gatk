package org.broadinstitute.hellbender.tools.walkers.vqsr;

import htsjdk.variant.vcf.VCFConstants;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

public class ApplyVQSRUnitTest extends GATKBaseTest {
    @Test
    public final void testGenerateFilterString() {
        final ApplyVQSR ar = new ApplyVQSR();
        ar.VQSLOD_CUTOFF = 0.0;
        Assert.assertTrue(ar.generateFilterString(5.0).equals(VCFConstants.PASSES_FILTERS_v4));
        Assert.assertTrue(ar.generateFilterString(-5.0).equals(ApplyVQSR.LOW_VQSLOD_FILTER_NAME));
    }
}
