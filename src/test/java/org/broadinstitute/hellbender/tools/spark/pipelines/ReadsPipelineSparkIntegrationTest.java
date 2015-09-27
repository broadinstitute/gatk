package org.broadinstitute.hellbender.tools.spark.pipelines;


import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.dataflow.datasources.RefAPISource;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.testng.Assert;
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
        final String localTestDir = publicTestDir + "org/broadinstitute/hellbender/tools/BQSR/";

        return new Object[][] {
                { localTestDir + "CEUTrio.HiSeq.WGS.b37.ch20.1m-1m1k.NA12878.bam",
                        RefAPISource.HS37D5_REF_ID, localTestDir + "dbsnp_132.b37.excluding_sites_after_129.chr17_69k_70k.vcf",
                        new File(localTestDir + "CEUTrio.HiSeq.WGS.b37.ch20.1m-1m1k.NA12878.bam") }
        };
    }

    @Test(dataProvider = "EndToEndTestData")
    public void testPipelineEndToEnd( final String inputBam, final String reference, final String knownSites, final File expectedOutput ) throws IOException, InterruptedException {
        File outFile = new File("output.bam");
        IOUtils.deleteRecursivelyOnExit(outFile);
        String output = "file:///" + outFile.getCanonicalPath();

        List<String> argv = new ArrayList<>();

        argv.addAll(Arrays.asList("-" + StandardArgumentDefinitions.INPUT_SHORT_NAME, inputBam,
                "-" + StandardArgumentDefinitions.REFERENCE_SHORT_NAME, RefAPISource.URL_PREFIX + reference,
                "-BQSRKnownVariants", knownSites,
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, output,
                "--apiKey", getDataflowTestApiKey()));
        int returnCode = runLocalSpark(argv);
        Assert.assertEquals(returnCode, 0);
    }
}