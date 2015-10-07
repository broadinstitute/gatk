package org.broadinstitute.hellbender.tools.dataflow.pipelines;


import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.dataflow.DataflowCommandLineProgramTest;
import org.broadinstitute.hellbender.engine.datasources.ReferenceAPISource;
import org.broadinstitute.hellbender.tools.IntegrationTestSpec;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class ReadsPreprocessingPipelineIntegrationTest extends DataflowCommandLineProgramTest {

    @DataProvider(name = "EndToEndTestData")
    public Object[][] getEndToEndTestData() {
        final String testDir = publicTestDir + "org/broadinstitute/hellbender/engine/";

        return new Object[][] {
                // This test case checks that the input and output bams are equal (which should be the case with stub tool implementations)
                { testDir + "reads_data_source_test1.bam", ReferenceAPISource.HS37D5_REF_ID, testDir + "feature_data_source_test.vcf", new File(testDir + "reads_data_source_test1.bam") }
        };
    }

    @Test(dataProvider = "EndToEndTestData", groups = {"cloud_todo"})
    public void testPipelineEndToEnd( final String inputBam, final String reference, final String knownSites, final File expectedOutput ) throws IOException {
        final File output = createTempFile("testPipelineEndToEnd_output", ".bam");
        List<String> argv = new ArrayList<>();

        argv.addAll(Arrays.asList("-" + StandardArgumentDefinitions.INPUT_SHORT_NAME, inputBam,
                                  "-" + StandardArgumentDefinitions.REFERENCE_SHORT_NAME, ReferenceAPISource.URL_PREFIX + reference,
                                  "-BQSRKnownVariants", knownSites,
                                  "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, output.getAbsolutePath()));
        argv.addAll(getStandardGCPArgumentsFromEnvironment());
        addDataflowRunnerArgs(argv);

        // Note: could use IntegrationTestSpec if we expand it a bit
        runCommandLine(argv);

        IntegrationTestSpec.assertEqualBamFiles(output, expectedOutput);
    }
}
