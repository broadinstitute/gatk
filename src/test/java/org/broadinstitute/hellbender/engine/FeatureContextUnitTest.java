package org.broadinstitute.hellbender.engine;

import htsjdk.tribble.Feature;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;

public class FeatureContextUnitTest extends BaseTest {

    @CommandLineProgramProperties(usage = "", usageShort = "")
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

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testNullIntervalWithNonNullFeatureManager() {
        // If we provide a backing data source, we must also provide a query interval
        try ( FeatureManager featureManager = new FeatureManager(new ArtificialFeatureContainingCommandLineProgram()) ) {
            FeatureContext featureContext = new FeatureContext(featureManager, null);
        }
    }

    @DataProvider(name = "EmptyFeatureContextDataProvider")
    public Object[][] getEmptyFeatureContextData() {
        // Default-constructed FeatureContexts and FeatureContexts constructed from null FeatureManagers
        // should behave as empty context objects.
        return new Object[][] {
                { new FeatureContext() },
                { new FeatureContext(null, null) },
                { new FeatureContext(null, hg19GenomeLocParser.createGenomeLoc("1", 1, 1) ) }
        };
    }

    @Test(dataProvider = "EmptyFeatureContextDataProvider")
    public void testFeatureContextWithNoBackingDataSource( final FeatureContext featureContext) {
        ArtificialFeatureContainingCommandLineProgram toolInstance = new ArtificialFeatureContainingCommandLineProgram();

        Assert.assertFalse(featureContext.hasBackingDataSource(), "Empty FeatureContext reports having a backing data source");
        Assert.assertTrue(featureContext.getValues(toolInstance.featureArgument).isEmpty(), "Empty FeatureContext should have returned an empty List from getValues()");
        Assert.assertTrue(featureContext.getValues(toolInstance.featureArgument, 1).isEmpty(), "Empty FeatureContext should have returned an empty List from getValues()");
        Assert.assertTrue(featureContext.getValues(Arrays.<FeatureInput<Feature>>asList(toolInstance.featureArgument)).isEmpty(), "Empty FeatureContext should have returned an empty List from getValues()");
        Assert.assertTrue(featureContext.getValues(Arrays.<FeatureInput<Feature>>asList(toolInstance.featureArgument), 1).isEmpty(), "Empty FeatureContext should have returned an empty List from getValues()");
    }

}
