package org.broadinstitute.hellbender.tools.spark.sv.integration;

import org.apache.hadoop.fs.Path;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.io.IOUtils;
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

    private static final File expectedSAMfile = IOUtils.getPath(SVIntegrationTestDataProvider.EXPECTED_ALIGNED_CONTIGS).toFile();

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

        String getCommandLine() {
            return  " -I " + bamLoc +
                    " -O "                    + outputDir + "/assemblies.sam" +
                    " --aligner-index-image " + alignerRefIndexImgLoc +
                    " --kmers-to-ignore " + kmerIgnoreListLoc +
                    " --breakpoint-intervals " + outputDir + "/intervals" +
                    " --fastq-dir "            + outputDir + "/fastq" +
                    " --target-link-file "      + outputDir + "/targetLinks.bedpe";
        }

        @Override
        public String toString() {
            return "FindBreakpointEvidenceSparkIntegrationTestArgs{" +
                    "bam-loc='" + bamLoc + '\'' +
                    ", kmer-ignore-list-loc='" + kmerIgnoreListLoc + '\'' +
                    ", aligner-ref-index-img-loc='" + alignerRefIndexImgLoc + '\'' +
                    ", output-dir='" + outputDir + '\'' +
                    '}';
        }
    }

    @DataProvider(name = "findBreakpointEvidenceSparkIntegrationTest")
    public Object[][] createTestData() throws IOException {

        List<Object[]> tests = new ArrayList<>();
        final File tempDirNew = BaseTest.createTempDir("forNew");
        tempDirNew.deleteOnExit();
        Files.createDirectories(Paths.get(tempDirNew.getAbsolutePath()+"/fastq"));
        tests.add(new Object[]{new FindBreakpointEvidenceSparkIntegrationTestArgs(SVIntegrationTestDataProvider.TEST_BAM,
                SVIntegrationTestDataProvider.KMER_KILL_LIST, SVIntegrationTestDataProvider.ALIGNER_INDEX_IMG,
                tempDirNew.getAbsolutePath())});

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "findBreakpointEvidenceSparkIntegrationTest", groups = "sv")
    public void testFindBreakpointRunnableLocal(final FindBreakpointEvidenceSparkIntegrationTestArgs params) throws IOException {

        final List<String> args = Arrays.asList( new ArgumentsBuilder().add(params.getCommandLine()).getArgsArray() );
        runCommandLine(args);

        final File actualSAMfile = new File(params.outputDir, "assemblies.sam");
        IntegrationTestSpec.assertEqualTextFiles(actualSAMfile, expectedSAMfile);
    }

    @Test(dataProvider = "findBreakpointEvidenceSparkIntegrationTest", groups = "sv")
    public void testFindBreakpointRunnableMiniCluster(final FindBreakpointEvidenceSparkIntegrationTestArgs params) throws Exception {

        MiniClusterUtils.runOnIsolatedMiniCluster(cluster -> {

            final List<String> argsToBeModified = Arrays.asList( new ArgumentsBuilder().add(params.getCommandLine()).getArgsArray() );
            final Path workingDirectory = MiniClusterUtils.getWorkingDir(cluster);

            int idx = 0;

            idx = argsToBeModified.indexOf("-I");
            Path path = new Path(workingDirectory, "hdfs.bam");
            File file = new File(argsToBeModified.get(idx+1));
            cluster.getFileSystem().copyFromLocalFile(new Path(file.toURI()), path);
            argsToBeModified.set(idx+1, path.toUri().toString());

            idx = argsToBeModified.indexOf("--kmers-to-ignore");
            path = new Path(workingDirectory, "dummy.kill.kmers");
            file = new File(argsToBeModified.get(idx+1));
            cluster.getFileSystem().copyFromLocalFile(new Path(file.toURI()), path);
            argsToBeModified.set(idx+1, path.toUri().toString());

            // outputs, prefix with hdfs address
            idx = argsToBeModified.indexOf("-O");
            path = new Path(workingDirectory, "assemblies.sam");
            argsToBeModified.set(idx+1, path.toUri().toString());

            idx = argsToBeModified.indexOf("--breakpoint-intervals");
            path = new Path(workingDirectory, "intervals");
            argsToBeModified.set(idx+1, path.toUri().toString());

            idx = argsToBeModified.indexOf("--fastq-dir");
            path = new Path(workingDirectory, "fastq");
            argsToBeModified.set(idx+1, path.toUri().toString());

            runCommandLine(argsToBeModified);

            final File actualSAMfile = BaseTest.createTempFile("assemblies", ".sam");
            cluster.getFileSystem().copyToLocalFile(new Path(workingDirectory, "assemblies.sam"), new Path(actualSAMfile.toURI()));
            IntegrationTestSpec.assertEqualTextFiles(actualSAMfile, expectedSAMfile);
        });
    }
}
