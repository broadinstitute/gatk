package org.broadinstitute.hellbender.tools.spark.sv.integration;

import org.apache.hadoop.fs.Path;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.test.MiniClusterUtils;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class DiscoverVariantsFromContigAlignmentsSAMSparkIntegrationTest extends CommandLineProgramTest {

    private static final class DiscoverVariantsFromContigAlignmentsSAMSparkIntegrationTestArgs {
        final String outputDir;

        DiscoverVariantsFromContigAlignmentsSAMSparkIntegrationTestArgs(final String outputDir){
            this.outputDir = outputDir;
        }

        String getCommandLineNoApiKey() {
            return  " -R " + SVIntegrationTestDataProvider.reference_2bit +
                    " -I " + SVIntegrationTestDataProvider.TEST_CONTIG_SAM +
                    " -O " + outputDir + "/variants.vcf" +
                    " --fastaReference " + SVIntegrationTestDataProvider.reference;
        }

    }

    @DataProvider(name = "discoverVariantsFromContigAlignmentsSparkIntegrationTest")
    public Object[][] createTestData() throws IOException {
        List<Object[]> tests = new ArrayList<>();
        final File tempDirLeft = BaseTest.createTempDir("forLeft");
        tempDirLeft.deleteOnExit();
        tests.add(new Object[]{new DiscoverVariantsFromContigAlignmentsSAMSparkIntegrationTestArgs(tempDirLeft.getAbsolutePath())});
        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "discoverVariantsFromContigAlignmentsSparkIntegrationTest", groups = "sv")
    public void testDiscoverVariantsRunnableLocal(final DiscoverVariantsFromContigAlignmentsSAMSparkIntegrationTestArgs params) throws Exception {

        final List<String> args = Arrays.asList( new ArgumentsBuilder().add(params.getCommandLineNoApiKey()).getArgsArray() );
        runCommandLine(args);
        StructuralVariationDiscoveryPipelineSparkIntegrationTest.svDiscoveryVCFEquivalenceTest(args.get(args.indexOf("-O")+1), SVIntegrationTestDataProvider.EXPECTED_SIMPLE_DEL_VCF, Arrays.asList("ALIGN_LENGTHS", "CTG_NAMES"), false);
    }

    @Test(dataProvider = "discoverVariantsFromContigAlignmentsSparkIntegrationTest", groups = "sv")
    public void testDiscoverVariantsRunnableMiniCluster(final DiscoverVariantsFromContigAlignmentsSAMSparkIntegrationTestArgs params) throws Exception {

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
            final String vcfOnHDFS = path.toUri().toString();
            argsToBeModified.set(idx+1, vcfOnHDFS);

            runCommandLine(argsToBeModified);
            StructuralVariationDiscoveryPipelineSparkIntegrationTest.svDiscoveryVCFEquivalenceTest(vcfOnHDFS, SVIntegrationTestDataProvider.EXPECTED_SIMPLE_DEL_VCF, Arrays.asList("ALIGN_LENGTHS", "CTG_NAMES"), true);
        });
    }
}
