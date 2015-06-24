package org.broadinstitute.hellbender.tools.walkers.bqsr;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.IntegrationTestSpec;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public final class ApplyWholeBQSRDataflowIntegrationTest extends CommandLineProgramTest {
    private static class ABQSRTest {
        final String bam;
        final String args;
        final String expectedFile;
        final String reference;
        final String knownSites;

        private ABQSRTest(String reference, String bam, String knownSites, String args, String expectedFile) {
            this.reference = reference;
            this.bam = bam;
            this.knownSites = knownSites;
            this.args = args;
            this.expectedFile = expectedFile;
        }

        @Override
        public String toString() {
            return String.format("ApplyBQSR(args='%s')", args);
        }
    }

    final String resourceDir = getTestDataDir() + "/" + "BQSR" + "/";
    final String hiSeqBam = resourceDir + "HiSeq.1mb.1RG.2k_lines.alternate_allaligned.bam";
    private String getResourceDir(){
        return resourceDir;
    }

    @DataProvider(name = "ApplyBQSRTest")
    public Object[][] createABQSRTestData() {
        final String hg18Reference = publicTestDir + "human_g1k_v37.chr17_1Mb.fasta";
        final String bam = getResourceDir() + "na.bam";
        final String dbSNPb37 =  getResourceDir() + "dbsnp_132.b37.excluding_sites_after_129.chr17_69k_70k.vcf";


        List<Object[]> tests = new ArrayList<>();

        tests.add(new Object[]{new ABQSRTest(hg18Reference, bam, dbSNPb37, "", resourceDir + "expected.na.bam")});

        // TODO: add test inputs with some unaligned reads

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "ApplyBQSRTest")
    public void testBQSR(ABQSRTest params) throws IOException {
        String dummy = BaseTest.createTempFile("temp-recaltables","tmp").getAbsolutePath();
        String args =
                " -R " + params.reference +
                " -I " + params.bam +
                (params.knownSites.isEmpty() ? "": " -knownSites " + params.knownSites) +
                // this argument is required, even though we don't really want it
                " --RECAL_TABLE_FILE " + dummy +
                params.args +
                " -O %s";
        ArgumentsBuilder ab = new ArgumentsBuilder().add(args);
        addDataflowRunnerArgs(ab);
        IntegrationTestSpec spec = new IntegrationTestSpec(
                ab.getString(),
                Arrays.asList(params.expectedFile));
        spec.executeTest("testBQSR-" + params.args, this);
    }

}
