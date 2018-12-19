package org.broadinstitute.hellbender.engine;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.cmdline.TestProgramGroup;
import org.testng.Assert;
import org.testng.annotations.Test;

public class MultiplePassVariantWalkerUnitTest extends GATKBaseTest {
    @CommandLineProgramProperties(
            summary = "An example subclass of MultiplePassVariantWalker",
            oneLineSummary = "An example subclass of MultiplePassVariantWalker",
            programGroup = TestProgramGroup.class,
            omitFromCommandLine = true
    )
    private static class DummyMultiplePassVariantWalker extends TwoPassVariantWalker {
        public int firstPass = 0;
        public int secondPass = 0;
        public boolean visitedAfterFirstPass = false;

        @Override
        protected void firstPassApply(VariantContext variant, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext) {
            firstPass++;
        }

        @Override
        protected void afterFirstPass(){
            visitedAfterFirstPass = true;
        }

        @Override
        protected void secondPassApply(VariantContext variant, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext) {
            secondPass++;
        }
    }

    @Test
    public void testTwoPassTraversal() {
        final DummyMultiplePassVariantWalker walker = new DummyMultiplePassVariantWalker();
        final String testVcf = "src/test/resources/org/broadinstitute/hellbender/tools/walkers/variantutils/VariantsToTable/multiallelic.vcf";

        final String[] args = { "-V", testVcf };

        walker.instanceMain(args);

        final int expectedNumberOfVariantContexts = 52;
        Assert.assertEquals(walker.firstPass, expectedNumberOfVariantContexts);
        Assert.assertEquals(walker.secondPass, expectedNumberOfVariantContexts);
        Assert.assertTrue(walker.visitedAfterFirstPass);
    }

}