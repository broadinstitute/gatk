package org.broadinstitute.hellbender.tools.spark.sv;

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

public class StructuralVariationDiscoveryPipelineSparkIntegrationTest extends CommandLineProgramTest {

    private static final class StructuralVariationDiscoveryPipelineSparkIntegrationTestArgs {
        final String bamLoc;
        final String kmerIgnoreListLoc;
        final String alignerRefIndexImgLoc;
        final String outputDir;


        StructuralVariationDiscoveryPipelineSparkIntegrationTestArgs(final String bamLoc,
                                                                     final String kmerIgnoreListLoc,
                                                                     final String alignerRefIndexImgLoc,
                                                                     final String outputDir) {
            this.bamLoc = bamLoc;
            this.kmerIgnoreListLoc = kmerIgnoreListLoc;
            this.alignerRefIndexImgLoc = alignerRefIndexImgLoc;
            this.outputDir = outputDir;
        }

        String getCommandLineNoApiKey() {
            return  " -R " + SVIntegrationTestDataProvider.reference_2bit +
                    " -I " + bamLoc +
                    " -O " + outputDir        + "/variants.vcf" +
                    " --fastaReference " + SVIntegrationTestDataProvider.reference +
                    " --alignerIndexImage " + alignerRefIndexImgLoc +
                    " --kmersToIgnore " + kmerIgnoreListLoc +
                    " --contigSAMFile "       + outputDir + "/assemblies.sam" +
                    " --breakpointIntervals " + outputDir + "/intervals" +
                    " --fastqDir "            + outputDir + "/fastq";
        }

        String getCommandLine() {
            return  getCommandLineNoApiKey() +
                    " --apiKey " + getGCPTestApiKey();
        }
    }

    @DataProvider(name = "svDiscoverPipelineSparkIntegrationTest")
    public Object[][] createTestData() throws IOException {
        List<Object[]> tests = new ArrayList<>();
        final File tempDirLeft = BaseTest.createTempDir("forLeft");
        tempDirLeft.deleteOnExit();
        Files.createDirectories(Paths.get(tempDirLeft.getAbsolutePath()+"/fastq"));
        tests.add(new Object[]{new StructuralVariationDiscoveryPipelineSparkIntegrationTest.StructuralVariationDiscoveryPipelineSparkIntegrationTestArgs(SVIntegrationTestDataProvider.TEST_BAM_LEFT, SVIntegrationTestDataProvider.KMER_KILL_LIST, SVIntegrationTestDataProvider.ALIGNER_INDEX_IMG, tempDirLeft.getAbsolutePath())});
        final File tempDirRight = BaseTest.createTempDir("forRight");
        tempDirRight.deleteOnExit();
        Files.createDirectories(Paths.get(tempDirRight.getAbsolutePath()+"/fastq"));
        tests.add(new Object[]{new StructuralVariationDiscoveryPipelineSparkIntegrationTest.StructuralVariationDiscoveryPipelineSparkIntegrationTestArgs(SVIntegrationTestDataProvider.TEST_BAM_RIGHT, SVIntegrationTestDataProvider.KMER_KILL_LIST, SVIntegrationTestDataProvider.ALIGNER_INDEX_IMG, tempDirRight.getAbsolutePath())});

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "svDiscoverPipelineSparkIntegrationTest", groups = "sv")
    public void testSVDiscoverPipelineRunnableLocal(final StructuralVariationDiscoveryPipelineSparkIntegrationTest.StructuralVariationDiscoveryPipelineSparkIntegrationTestArgs params) throws IOException {
        new IntegrationTestSpec(
                new ArgumentsBuilder().add(params.getCommandLineNoApiKey()).getString(),
                SVIntegrationTestDataProvider.dummyExpectedFileNames)
                .executeTest("testSVDiscoverPipelineRunnableLocal-", this);
    }

    @Test(dataProvider = "svDiscoverPipelineSparkIntegrationTest", groups = "sv")
    public void testSVDiscoverPipelineRunnableMiniCluster(final StructuralVariationDiscoveryPipelineSparkIntegrationTest.StructuralVariationDiscoveryPipelineSparkIntegrationTestArgs params) throws Exception {

        MiniClusterUtils.runOnIsolatedMiniCluster(cluster -> {

            final List<String> argsToBeModified = Arrays.asList( new ArgumentsBuilder().add(params.getCommandLineNoApiKey()).getArgsArray() );
            final Path workingDirectory = MiniClusterUtils.getWorkingDir(cluster);

            int idx = 0;

            // inputs, copy to mini cluster
            idx = argsToBeModified.indexOf("-I");
            Path path = new Path(workingDirectory, "hdfs.bam");
            File file = new File(argsToBeModified.get(idx+1));
            cluster.getFileSystem().copyFromLocalFile(new Path(file.toURI()), path);
            argsToBeModified.set(idx+1, path.toUri().toString());

            idx = argsToBeModified.indexOf("-R");
            path = new Path(workingDirectory, "reference.2bit");
            file = new File(argsToBeModified.get(idx+1));
            cluster.getFileSystem().copyFromLocalFile(new Path(file.toURI()), path);
            argsToBeModified.set(idx+1, path.toUri().toString());

            idx = argsToBeModified.indexOf("--fastaReference");
            path = new Path(workingDirectory, "reference.fasta");
            file = new File(argsToBeModified.get(idx+1));
            cluster.getFileSystem().copyFromLocalFile(new Path(file.toURI()), path);
            argsToBeModified.set(idx+1, path.toUri().toString());

            idx = argsToBeModified.indexOf("--kmersToIgnore");
            path = new Path(workingDirectory, "dummy.kill.kmers");
            file = new File(argsToBeModified.get(idx+1));
            cluster.getFileSystem().copyFromLocalFile(new Path(file.toURI()), path);
            argsToBeModified.set(idx+1, path.toUri().toString());

            path = new Path(workingDirectory, "reference.fasta.fai");
            cluster.getFileSystem().copyFromLocalFile(new Path(SVIntegrationTestDataProvider.reference_fai.toURI()), path);
            path = new Path(workingDirectory, "reference.dict");
            cluster.getFileSystem().copyFromLocalFile(new Path(SVIntegrationTestDataProvider.reference_dict.toURI()), path);

            // outputs, prefix with hdfs address
            idx = argsToBeModified.indexOf("-O");
            path = new Path(workingDirectory, "variants.vcf");
            argsToBeModified.set(idx+1, path.toUri().toString());

            idx = argsToBeModified.indexOf("--contigSAMFile");
            path = new Path(workingDirectory, "assemblies.sam");
            argsToBeModified.set(idx+1, path.toUri().toString());

            idx = argsToBeModified.indexOf("--breakpointIntervals");
            path = new Path(workingDirectory, "intervals");
            argsToBeModified.set(idx+1, path.toUri().toString());

            idx = argsToBeModified.indexOf("--fastqDir");
            path = new Path(workingDirectory, "fastq");
            argsToBeModified.set(idx+1, path.toUri().toString());

            new IntegrationTestSpec(String.join(" ", argsToBeModified), SVIntegrationTestDataProvider.dummyExpectedFileNames)
                    .executeTest("testSVDiscoverPipelineRunnableMiniCluster-", this);
        });
    }
}
