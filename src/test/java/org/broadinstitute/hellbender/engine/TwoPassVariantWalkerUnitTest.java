package org.broadinstitute.hellbender.engine;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.TestProgramGroup;
import org.testng.Assert;
import org.testng.annotations.Test;

public class TwoPassVariantWalkerUnitTest {
    @CommandLineProgramProperties(
            summary = "An example subclass of TwoPassVariantWalker",
            oneLineSummary = "An example subclass of TwoPassVariantWalker",
            programGroup = TestProgramGroup.class
    )
    private static class DummyTwoPassVariantWalker extends TwoPassVariantWalker {
        public int firstPass = 0;
        public int secondPass = 0;
        public boolean visitedAfterFirstPass = false;

        @Override
        protected void firstPassApply(VariantContext variant, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext) {
            firstPass++;
        }

        @Override
        protected void secondPassApply(VariantContext variant, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext) {
            secondPass++;
        }

        @Override
        protected void afterFirstPass(){
            visitedAfterFirstPass = true;
        }
    }

    @Test
    public void testNum() {
        final DummyTwoPassVariantWalker walker = new DummyTwoPassVariantWalker();
        final String testVcf = "src/test/resources/org/broadinstitute/hellbender/tools/walkers/variantutils/VariantsToTable/multiallelic.vcf";

        final String[] args = { "-V", testVcf };

        walker.instanceMain(args);

        final int expectedNumberOfVariantContexts = 52;
        Assert.assertEquals(walker.firstPass, expectedNumberOfVariantContexts);
        Assert.assertEquals(walker.secondPass, expectedNumberOfVariantContexts);
        Assert.assertTrue(walker.visitedAfterFirstPass);
    }

}