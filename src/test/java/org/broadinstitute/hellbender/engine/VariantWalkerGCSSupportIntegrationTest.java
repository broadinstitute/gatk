package org.broadinstitute.hellbender.engine;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.tools.walkers.variantutils.SelectVariants;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.IOException;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * Tests to prove that we can access and query inputs on Google Cloud Storage (GCS) in VariantWalkers.
 */
public class VariantWalkerGCSSupportIntegrationTest extends CommandLineProgramTest {

    private static final String TEST_VCF_ON_GCS = "org/broadinstitute/hellbender/engine/dbsnp_138.b37.20.10000000-10010000.vcf";
    private static final String TEST_BGZIPPED_VCF_ON_GCS = "org/broadinstitute/hellbender/engine/dbsnp_138.b37.20.10000000-10010000.vcf.block.gz";
    private static final String EXPECTED_OUTPUT_DIR = publicTestDir + "org/broadinstitute/hellbender/engine/GCSTests/";

    @Override
    public String getTestedToolName() {
        return SelectVariants.class.getSimpleName();
    }

    @DataProvider(name = "GCSTestCases")
    public Object[][] getGCSTestCases() {
        final SimpleInterval singleInterval = new SimpleInterval("20", 10004000, 10006000);
        final List<SimpleInterval> multipleIntervals = Arrays.asList(singleInterval, new SimpleInterval("20", 10008000, 10009000));

        final String EXPECTED_WHOLE_FILE_RESULTS = EXPECTED_OUTPUT_DIR + "expected_VariantWalkerGCSSupportIntegrationTest_vcf_wholefile.vcf";
        final String EXPECTED_SINGLE_INTERVAL_RESULTS = EXPECTED_OUTPUT_DIR + "expected_VariantWalkerGCSSupportIntegrationTest_vcf_single_interval.vcf";
        final String EXPECTED_MULTIPLE_INTERVALS_RESULTS = EXPECTED_OUTPUT_DIR + "expected_VariantWalkerGCSSupportIntegrationTest_vcf_multiple_intervals.vcf";

        return new Object[][] {
                { TEST_VCF_ON_GCS, null, EXPECTED_WHOLE_FILE_RESULTS },
                { TEST_BGZIPPED_VCF_ON_GCS, null, EXPECTED_WHOLE_FILE_RESULTS },
                { TEST_VCF_ON_GCS, Collections.singletonList(singleInterval), EXPECTED_SINGLE_INTERVAL_RESULTS },
                { TEST_BGZIPPED_VCF_ON_GCS, Collections.singletonList(singleInterval), EXPECTED_SINGLE_INTERVAL_RESULTS },
                { TEST_VCF_ON_GCS, multipleIntervals, EXPECTED_MULTIPLE_INTERVALS_RESULTS },
                { TEST_BGZIPPED_VCF_ON_GCS, multipleIntervals, EXPECTED_MULTIPLE_INTERVALS_RESULTS }
        };
    }

    @Test(dataProvider = "GCSTestCases", groups = {"bucket"})
    public void testReadVCFOnGCS( final String vcf, final List<SimpleInterval> intervals, final String expectedOutput ) throws IOException {
        String intervalArg = "";
        if ( intervals != null ) {
            final StringBuilder intervalArgBuilder = new StringBuilder();
            for ( final SimpleInterval interval : intervals ) {
                intervalArgBuilder.append(" -L ");
                intervalArgBuilder.append(interval.toString());
            }
            intervalArg = intervalArgBuilder.toString();
        }

        final IntegrationTestSpec testSpec = new IntegrationTestSpec(
                " -V " + getGCPTestInputPath() + vcf +
                intervalArg +
                " -O %s "+ " --" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE + " false",
                Collections.singletonList(expectedOutput)
        );
        testSpec.executeTest("testReadVCFOnGCS", this);
    }
}
