package org.broadinstitute.hellbender.engine;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;

public class VariantWalkerBaseUnitTest extends CommandLineProgramTest {
    final String baseVariants = packageRootTestDir + "engine/feature_data_source_test.vcf";

    private static class TestVariantWalker extends VariantWalker {
        public TestVariantWalker(String drivingVariants) {
            this.drivingVariantFile = drivingVariants;
        }

        @Override
        public void apply(VariantContext variant, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext) {

        }
    }

    @Test
    public void testGetIntervals() {

        final TestVariantWalker tool = new TestVariantWalker(baseVariants);
        tool.instanceMain(new String[]{
                "-V", baseVariants,
                "-R", hg19MiniReference,
                "-L", "1:21-21"
        });

        Assert.assertEquals(tool.getIntervals().size(), 1);

        final TestVariantWalker tool2 = new TestVariantWalker(baseVariants);
        tool2.instanceMain(new String[]{
                "-V", baseVariants,
                "-R", hg19MiniReference
        });

        Assert.assertEquals(tool2.getIntervals().size(), 4);
    }
}
