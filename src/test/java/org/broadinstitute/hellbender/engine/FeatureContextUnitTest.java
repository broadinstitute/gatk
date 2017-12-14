package org.broadinstitute.hellbender.engine;

import org.broadinstitute.hellbender.utils.SimpleInterval;
import htsjdk.tribble.Feature;
import htsjdk.variant.vcf.VCFHeader;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.cmdline.TestProgramGroup;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;

public final class FeatureContextUnitTest extends GATKBaseTest {

    @CommandLineProgramProperties(summary = "", oneLineSummary = "", programGroup = TestProgramGroup.class)
    private static class ArtificialFeatureContainingCommandLineProgram extends CommandLineProgram {
        @Argument(fullName = "featureArgument", shortName = "f")
        FeatureInput<Feature> featureArgument;

        public ArtificialFeatureContainingCommandLineProgram() {
            featureArgument = new FeatureInput<>(publicTestDir + "org/broadinstitute/hellbender/engine/feature_data_source_test.vcf");
        }

        @Override
        protected Object doWork() {
            return null;
        }
    }

    @DataProvider(name = "EmptyFeatureContextDataProvider")
    public Object[][] getEmptyFeatureContextData() {
        // Default-constructed FeatureContexts and FeatureContexts constructed from null FeatureManagers
        // or intervals should behave as empty context objects.
        return new Object[][] {
                { new FeatureContext() },
                { new FeatureContext(null, null) },
                { new FeatureContext(null, new SimpleInterval("1", 1, 1) ) },
                { new FeatureContext(new FeatureManager(new ArtificialFeatureContainingCommandLineProgram()), null) }
        };
    }

    @Test(dataProvider = "EmptyFeatureContextDataProvider")
    public void testEmptyFeatureContext( final FeatureContext featureContext) {
        ArtificialFeatureContainingCommandLineProgram toolInstance = new ArtificialFeatureContainingCommandLineProgram();

        Assert.assertFalse(featureContext.hasBackingDataSource() && featureContext.getInterval() != null,
                           "Empty FeatureContext reports having both a backing data source and an interval");
        Assert.assertTrue(featureContext.getValues(toolInstance.featureArgument).isEmpty(), "Empty FeatureContext should have returned an empty List from getValues()");
        Assert.assertTrue(featureContext.getValues(toolInstance.featureArgument, 1).isEmpty(), "Empty FeatureContext should have returned an empty List from getValues()");
        Assert.assertTrue(featureContext.getValues(Arrays.<FeatureInput<Feature>>asList(toolInstance.featureArgument)).isEmpty(), "Empty FeatureContext should have returned an empty List from getValues()");
        Assert.assertTrue(featureContext.getValues(Arrays.<FeatureInput<Feature>>asList(toolInstance.featureArgument), 1).isEmpty(), "Empty FeatureContext should have returned an empty List from getValues()");
    }

    @Test
    public void testGetHeader() {
        final ArtificialFeatureContainingCommandLineProgram toolInstance = new ArtificialFeatureContainingCommandLineProgram();
        try (final FeatureManager featureManager = new FeatureManager(toolInstance)) {
            final FeatureContext featureContext = new FeatureContext(featureManager, new SimpleInterval("1", 1, 1));
            final Object header = featureContext.getHeader(toolInstance.featureArgument);
            Assert.assertTrue(header instanceof VCFHeader, "Header for " + toolInstance.featureArgument.getFeaturePath() +
                    " not a VCFHeader");
        }
    }
}
