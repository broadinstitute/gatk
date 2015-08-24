package org.broadinstitute.hellbender.tools.spark.pipelines;


import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.dataflow.datasources.RefAPISource;
import org.broadinstitute.hellbender.tools.IntegrationTestSpec;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class ReadsPipelineSparkIntegrationTest extends CommandLineProgramTest {

    @DataProvider(name = "EndToEndTestData")
    public Object[][] getEndToEndTestData() {
        final String testDir = "gs://davidada-hellbender-ctd/resources/";
        final String localTestDir = publicTestDir + "org/broadinstitute/hellbender/tools/BQSR/";

        return new Object[][] {
                // This test case checks that the input and output bams are equal (which should be the case with stub tool implementations)
                { testDir + "CEUTrio.HiSeq.WGS.b37.ch20.1m-1m1k.NA12878.bam",
                        RefAPISource.HS37D5_REF_ID, testDir + "dbsnp_132.b37.chr17_69k_70k.vcf",
                        new File(localTestDir + "CEUTrio.HiSeq.WGS.b37.ch20.1m-1m1k.NA12878.bam") }
        };
    }

    @Test(dataProvider = "EndToEndTestData", groups = {"cloud_todo"})
    public void testPipelineEndToEnd( final String inputBam, final String reference, final String knownSites, final File expectedOutput ) throws IOException, InterruptedException {
        String output = "gs://davidada-hellbender-ctd/output/output.bam"; //createTempFile("testPipelineEndToEnd_output", ".bam");
        List<String> argv = new ArrayList<>();

        argv.addAll(Arrays.asList("-" + StandardArgumentDefinitions.INPUT_SHORT_NAME, inputBam,
                                  "-" + StandardArgumentDefinitions.REFERENCE_SHORT_NAME, RefAPISource.URL_PREFIX + reference,
                                  "-BQSRKnownVariants", knownSites,
                                  "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, output,
                                  "--apiKey", getDataflowTestApiKey(),
                                  "--sparkMaster", "spark://broad-dsde-dev-m:7077"));
        //argv.addAll(getStandardDataflowArgumentsFromEnvironment());

        // 1) gradle shadowjar (if necessary)
        // 2) gcloud compute copy-files \
        // build/libs/hellbender_branch-all-GATK.4.alpha-*-SNAPSHOT-spark.jar broad-dsde-dev-m:~/ --zone "us-central1-a"
        // 3) ssh command to run (some key-gen happens, not sure how to handle that).
        // 4)
        runSparkCommandLine(argv);

        //IntegrationTestSpec.assertEqualBamFiles(output, expectedOutput);
    }
}
