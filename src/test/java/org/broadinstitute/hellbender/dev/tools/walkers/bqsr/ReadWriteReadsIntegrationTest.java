package org.broadinstitute.hellbender.dev.tools.walkers.bqsr;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.IntegrationTestSpec;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/*
 * Check that we can round-trip the files we are testing with.
 */
public final class ReadWriteReadsIntegrationTest extends CommandLineProgramTest {

    private final static String THIS_TEST_FOLDER = "org/broadinstitute/hellbender/tools/BQSR/";

    final String resourceDir = getTestDataDir() + "/" + "BQSR" + "/";
    final String hiSeqBam = resourceDir + "HiSeq.1mb.1RG.2k_lines.alternate_allaligned.bam";
    final String naBam = resourceDir + "NA12878.chr17_69k_70k.dictFix.bam";

    public static final class Bams {
        public final String localBam;
        public final String remoteBam;
        public Bams(String local, String remote) {
            this.localBam = local;
            this.remoteBam = remote;
        }
    }


    @DataProvider(name = "ApplyTest")
    public Object[][] createTestData() {
        List<Object[]> tests = new ArrayList<>();

        tests.add(new Object[]{hiSeqBam});

        // TODO: add test inputs with some unaligned reads

        return tests.toArray(new Object[][]{});
    }

    @DataProvider(name = "ApplyTest_GCS")
    public Object[][] createTestData_GCS() {
        final String resourceDirGCS = getDataflowTestInputPath() + THIS_TEST_FOLDER;
        final String hiSeqBamGCS = resourceDirGCS + "HiSeq.1mb.1RG.2k_lines.alternate_allaligned.bam";
        List<Object[]> tests = new ArrayList<>();

        // This shouldn't fail, but it does.
        tests.add(new Object[]{new Bams(hiSeqBam, hiSeqBamGCS)});

        return tests.toArray(new Object[][]{});
    }


    @Test(dataProvider = "ApplyTest")
    public void testRW(String bam) throws IOException {
        String args =
                " -I " + bam +
                " -O %s";
        ArgumentsBuilder ab = new ArgumentsBuilder().add(args);
        addDataflowRunnerArgs(ab);
        IntegrationTestSpec spec = new IntegrationTestSpec(
                ab.getString(),
                Arrays.asList(bam));
        spec.executeTest("testRW-" + bam, this);
    }

    @Test(dataProvider = "ApplyTest_GCS")
    public void testRW_GCS(Bams bams) throws IOException {
        String args =
                " -I " + bams.remoteBam +
                " -O %s";
        ArgumentsBuilder ab = new ArgumentsBuilder().add(args);
        addDataflowRunnerArgs(ab);
        IntegrationTestSpec spec = new IntegrationTestSpec(
                ab.getString(),
                Arrays.asList(bams.localBam));
        spec.executeTest("testRW_GCS-" + bams.localBam, this);
    }

}
