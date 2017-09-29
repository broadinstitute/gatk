package org.broadinstitute.hellbender.tools.spark.sv.integration;

import org.apache.hadoop.fs.Path;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.test.IntegrationTestSpec;
import org.broadinstitute.hellbender.utils.test.MiniClusterUtils;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Integration test on the SV pipeline as it exists right now [2017-03-06]
 */
public class FindBreakpointEvidenceSparkIntegrationTest extends CommandLineProgramTest {

    private static final class FindBreakpointEvidenceSparkIntegrationTestArgs {
        final String bamLoc;
        final String kmerIgnoreListLoc;
        final String alignerRefIndexImgLoc;
        final String outputDir;

        FindBreakpointEvidenceSparkIntegrationTestArgs(final String bamLoc,
                                                       final String kmerIgnoreListLoc,
                                                       final String alignerRefIndexImgLoc,
                                                       final String outputDir) {
            this.bamLoc = bamLoc;
            this.kmerIgnoreListLoc = kmerIgnoreListLoc;
            this.alignerRefIndexImgLoc = alignerRefIndexImgLoc;
            this.outputDir = outputDir;
        }

        String getCommandLineNoApiKey() {
            return  " -I " + bamLoc +
                    " -O "                    + outputDir + "/assemblies.sam" +
                    " --alignerIndexImage " + alignerRefIndexImgLoc +
                    " --kmersToIgnore " + kmerIgnoreListLoc +
                    " --breakpointIntervals " + outputDir + "/intervals" +
                    " --fastqDir "            + outputDir + "/fastq";
        }
    }

    @DataProvider(name = "findBreakpointEvidenceSparkIntegrationTest")
    public Object[][] createTestData() throws IOException {

        List<Object[]> tests = new ArrayList<>();
        final File tempDirLeft = BaseTest.createTempDir("forLeft");
        tempDirLeft.deleteOnExit();
        Files.createDirectories(Paths.get(tempDirLeft.getAbsolutePath()+"/fastq"));
        tests.add(new Object[]{new FindBreakpointEvidenceSparkIntegrationTestArgs(SVIntegrationTestDataProvider.TEST_BAM_LEFT, SVIntegrationTestDataProvider.KMER_KILL_LIST, SVIntegrationTestDataProvider.ALIGNER_INDEX_IMG, tempDirLeft.getAbsolutePath())});
        final File tempDirRight = BaseTest.createTempDir("forRight");
        tempDirRight.deleteOnExit();
        Files.createDirectories(Paths.get(tempDirRight.getAbsolutePath()+"/fastq"));
        tests.add(new Object[]{new FindBreakpointEvidenceSparkIntegrationTestArgs(SVIntegrationTestDataProvider.TEST_BAM_RIGHT, SVIntegrationTestDataProvider.KMER_KILL_LIST, SVIntegrationTestDataProvider.ALIGNER_INDEX_IMG, tempDirRight.getAbsolutePath())});

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "findBreakpointEvidenceSparkIntegrationTest", groups = "sv")
    public void testFindBreakpointRunnableLocal(final FindBreakpointEvidenceSparkIntegrationTestArgs params) throws IOException {

        new IntegrationTestSpec(
                new ArgumentsBuilder().add(params.getCommandLineNoApiKey()).getString(),
                SVIntegrationTestDataProvider.dummyExpectedFileNames)
                .executeTest("testFindBreakpointEvidenceSparkRunnableLocal-", this);
    }

    @Test(dataProvider = "findBreakpointEvidenceSparkIntegrationTest", groups = "sv")
    public void testFindBreakpointRunnableMiniCluster(final FindBreakpointEvidenceSparkIntegrationTestArgs params) throws Exception {

        MiniClusterUtils.runOnIsolatedMiniCluster(cluster -> {

            final List<String> argsToBeModified = Arrays.asList( new ArgumentsBuilder().add(params.getCommandLineNoApiKey()).getArgsArray() );
            final Path workingDirectory = MiniClusterUtils.getWorkingDir(cluster);

            int idx = 0;

            idx = argsToBeModified.indexOf("-I");
            Path path = new Path(workingDirectory, "hdfs.bam");
            File file = new File(argsToBeModified.get(idx+1));
            cluster.getFileSystem().copyFromLocalFile(new Path(file.toURI()), path);
            argsToBeModified.set(idx+1, path.toUri().toString());

            idx = argsToBeModified.indexOf("--kmersToIgnore");
            path = new Path(workingDirectory, "dummy.kill.kmers");
            file = new File(argsToBeModified.get(idx+1));
            cluster.getFileSystem().copyFromLocalFile(new Path(file.toURI()), path);
            argsToBeModified.set(idx+1, path.toUri().toString());

            // outputs, prefix with hdfs address
            idx = argsToBeModified.indexOf("-O");
            path = new Path(workingDirectory, "assemblies.sam");
            argsToBeModified.set(idx+1, path.toUri().toString());

            idx = argsToBeModified.indexOf("--breakpointIntervals");
            path = new Path(workingDirectory, "intervals");
            argsToBeModified.set(idx+1, path.toUri().toString());

            idx = argsToBeModified.indexOf("--fastqDir");
            path = new Path(workingDirectory, "fastq");
            argsToBeModified.set(idx+1, path.toUri().toString());

            new IntegrationTestSpec(String.join(" ", argsToBeModified), SVIntegrationTestDataProvider.dummyExpectedFileNames)
                    .executeTest("testFindBreakpointEvidenceSparkRunnableMiniCluster-", this);
        });
    }
}
