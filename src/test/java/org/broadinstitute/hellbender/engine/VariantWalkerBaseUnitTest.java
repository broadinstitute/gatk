package org.broadinstitute.hellbender.engine;

import com.sun.istack.Nullable;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.TestProgramGroup;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class VariantWalkerBaseUnitTest extends CommandLineProgramTest {
    final String baseVariants = packageRootTestDir + "engine/feature_data_source_test.vcf";

    @CommandLineProgramProperties(programGroup = TestProgramGroup.class, oneLineSummary = "VariantWalkerBase Test Walker", summary = "This is a test walker for VariantWalkerBase")
    private static class TestVariantWalker extends VariantWalker {
        public TestVariantWalker(String drivingVariants) {
            this.drivingVariantFile = drivingVariants;
        }

        @Override
        public void apply(VariantContext variant, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext) {

        }
    }

    @DataProvider(name = "TestGetIntervalsProvider")
    public Object[][] getTestGetIntervalsProvider() {
        return new Object[][]{
                {"1:21-21", 1},
                {null, 4}
        };
    }

    @Test(dataProvider = "TestGetIntervalsProvider")
    public void testGetIntervals(@Nullable String intervals, int expected) {
        List<String> args = new ArrayList<>(Arrays.asList(
                "-V", baseVariants,
                "-R", hg19MiniReference
        ));

        if (intervals != null) {
            args.add("-L");
            args.add(intervals);
        }

        final TestVariantWalker tool = new TestVariantWalker(baseVariants);
        tool.instanceMain(args.toArray(new String[args.size()]));

        Assert.assertEquals(tool.getIntervals().size(), expected);
    }
}
