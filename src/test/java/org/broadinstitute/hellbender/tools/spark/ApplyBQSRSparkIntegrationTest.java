package org.broadinstitute.hellbender.tools.spark;

import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.broadinstitute.hellbender.tools.walkers.bqsr.AbstractApplyBQSRIntegrationTest;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

@Test(groups = "spark")
public final class ApplyBQSRSparkIntegrationTest extends AbstractApplyBQSRIntegrationTest {
    private final static String THIS_TEST_FOLDER = toolsTestDir + "BQSR/";

    @Override
    protected String getDevNull() throws IOException {
        File file = File.createTempFile("devnull", "");
        file.deleteOnExit();
        return file.getAbsolutePath();
    }

    @DataProvider(name = "ApplyBQSRTestGCS")
    public Object[][] createABQSRTestDataGCS() {
        final String resourceDirGCS = getGCPTestInputPath() + THIS_TEST_FOLDER;
        final String hiSeqBamGCS = resourceDirGCS + "HiSeq.1mb.1RG.2k_lines.alternate_allaligned.bam";

        List<Object[]> tests = new ArrayList<>();
        tests.add(new Object[]{new ABQSRTest(hiSeqBamGCS, null, ".bam", null, resourceDir + "expected.HiSeq.1mb.1RG.2k_lines.alternate_allaligned.recalibrated.DIQ.bam")});

        // TODO: add test inputs with some unaligned reads

        return tests.toArray(new Object[][]{});
    }

    //TODO: This is disabled because we can't read a google bucket as a hadoop file system outside of the dataproc environment yet
    //Renable when we've figured out how to setup the google hadoop fs connector
    @Test(dataProvider = "ApplyBQSRTestGCS", groups = {"spark", "bucket"}, enabled = false)
    public void testPR_GCS(ABQSRTest params) throws IOException {
        String args =
                " -I " + params.bam +
                        " --" + StandardArgumentDefinitions.BQSR_TABLE_LONG_NAME + " " + resourceDir + "HiSeq.20mb.1RG.table.gz " +
                        params.args +
                        " -O %s";
        ArgumentsBuilder ab = new ArgumentsBuilder().add(args);
        IntegrationTestSpec spec = new IntegrationTestSpec(
                ab.getString(),
                Arrays.asList(params.expectedFile));
        spec.executeTest("testPrintReads-" + params.args, this);
    }

    @Test(dataProvider = "ApplyBQSRTestGCS", groups = {"spark", "bucket"}, enabled = false)
    public void testPR_Cloud(ABQSRTest params) throws IOException {
        String args =
                " -I " + params.bam +
                        " --" + StandardArgumentDefinitions.BQSR_TABLE_LONG_NAME + " " + getGCPTestInputPath() + THIS_TEST_FOLDER + "HiSeq.20mb.1RG.table.gz " +
                        params.args +
                        " -O %s";
        IntegrationTestSpec spec = new IntegrationTestSpec(
                args,
                Arrays.asList(params.expectedFile));
        spec.executeTest("testPrintReads-" + params.args, this);
    }
}
