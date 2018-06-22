package org.broadinstitute.hellbender.tools.spark.pipelines;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

public final class PrintVariantsSparkIntegrationTest extends CommandLineProgramTest {

    @Override
    public String getTestedClassName() {
        return PrintVariantsSpark.class.getSimpleName();
    }

    @DataProvider
    public Object[][] gcsTestingData() {
        return new Object[][]{
            {"org/broadinstitute/hellbender/engine/dbsnp_138.b37.20.10000000-10010000.vcf", ".vcf", false},
            {"org/broadinstitute/hellbender/engine/dbsnp_138.b37.20.10000000-10010000.vcf", ".vcf", true},
        };
    }

    /**
     * Test the Spark code locally, including GCS access.
     *
     * For this to work, the settings in src/main/resources/core-site.xml must be correct,
     * and the project name and credential file it points to must be present.
     */
    @Test(dataProvider = "gcsTestingData", groups = "bucket")
    public void testGCSInputsAndOutputs(final String gcsInput, final String outputExtension,
        final boolean outputToGCS) {
        final String gcsInputPath = getGCPTestInputPath() + gcsInput;
        final String outputPrefix = outputToGCS ? getGCPTestStaging() : "testGCSInputsAndOutputs";
        final String outputPath = BucketUtils.getTempFilePath(outputPrefix, outputExtension);

        final ArgumentsBuilder argBuilder = new ArgumentsBuilder();
        argBuilder.addArgument("variant", gcsInputPath)
            .addArgument("output", outputPath);
        runCommandLine(argBuilder);
    }
}
