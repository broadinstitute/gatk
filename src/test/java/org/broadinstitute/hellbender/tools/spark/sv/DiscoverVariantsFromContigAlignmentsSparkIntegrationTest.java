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
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class DiscoverVariantsFromContigAlignmentsSparkIntegrationTest extends CommandLineProgramTest {

    private static final class DiscoverVariantsFromContigAlignmentsSparkIntegrationTestArgs {
        final String outputDir;

        DiscoverVariantsFromContigAlignmentsSparkIntegrationTestArgs(final String outputDir){
            this.outputDir = outputDir;
        }

        String getCommandLineNoApiKey() {
            return  " -R " + SVIntegrationTestDataProvider.reference_2bit +
                    " -I " + SVIntegrationTestDataProvider.TEST_CONTIG_SAM +
                    " -O " + outputDir + "/variants.vcf" +
                    " --fastaReference " + SVIntegrationTestDataProvider.reference;
        }

        String getCommandLine() {
            return  getCommandLineNoApiKey() +
                    " --apiKey " + getGCPTestApiKey();
        }
    }

    @DataProvider(name = "discoverVariantsFromContigAlignmentsSparkIntegrationTest")
    public Object[][] createTestData() throws IOException {
        List<Object[]> tests = new ArrayList<>();
        final File tempDirLeft = BaseTest.createTempDir("forLeft");
        tempDirLeft.deleteOnExit();
        tests.add(new Object[]{new DiscoverVariantsFromContigAlignmentsSparkIntegrationTestArgs(tempDirLeft.getAbsolutePath())});
        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "discoverVariantsFromContigAlignmentsSparkIntegrationTest", groups = "sv")
    public void testDiscoverVariantsRunnableLocal(final DiscoverVariantsFromContigAlignmentsSparkIntegrationTest.DiscoverVariantsFromContigAlignmentsSparkIntegrationTestArgs params) throws IOException {
        new IntegrationTestSpec(
                new ArgumentsBuilder().add(params.getCommandLineNoApiKey()).getString(),
                SVIntegrationTestDataProvider.dummyExpectedFileNames)
                .executeTest("testDiscoverVariantsRunnableLocal-", this);
    }

    @Test(dataProvider = "discoverVariantsFromContigAlignmentsSparkIntegrationTest", groups = "sv")
    public void testDiscoverVariantsRunnableMiniCluster(final DiscoverVariantsFromContigAlignmentsSparkIntegrationTest.DiscoverVariantsFromContigAlignmentsSparkIntegrationTestArgs params) throws Exception {

        MiniClusterUtils.runOnIsolatedMiniCluster(cluster -> {

            final List<String> argsToBeModified = Arrays.asList( new ArgumentsBuilder().add(params.getCommandLineNoApiKey()).getArgsArray() );
            final Path workingDirectory = MiniClusterUtils.getWorkingDir(cluster);

            int idx = 0;

            idx = argsToBeModified.indexOf("-I");
            Path path = new Path(workingDirectory, "hdfs.sam");
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

            path = new Path(workingDirectory, "reference.fasta.fai");
            cluster.getFileSystem().copyFromLocalFile(new Path(SVIntegrationTestDataProvider.reference_fai.toURI()), path);
            path = new Path(workingDirectory, "reference.dict");
            cluster.getFileSystem().copyFromLocalFile(new Path(SVIntegrationTestDataProvider.reference_dict.toURI()), path);

            // outputs, prefix with hdfs address
            idx = argsToBeModified.indexOf("-O");
            path = new Path(workingDirectory, "variants.vcf");
            argsToBeModified.set(idx+1, path.toUri().toString());

            new IntegrationTestSpec(String.join(" ", argsToBeModified), SVIntegrationTestDataProvider.dummyExpectedFileNames)
                    .executeTest("testDiscoverVariantsRunnableMiniCluster-", this);
        });
    }
}
