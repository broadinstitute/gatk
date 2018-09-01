package org.broadinstitute.hellbender.tools.spark.sv.integration;

import org.apache.hadoop.fs.Path;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.testutils.MiniClusterUtils;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class DiscoverVariantsFromContigAlignmentsSAMSparkIntegrationTest extends CommandLineProgramTest {

    static final List<String> annotationsToIgnoreWhenComparingVariants =
            Arrays.asList(GATKSVVCFConstants.ALIGN_LENGTHS,
                    GATKSVVCFConstants.CONTIG_NAMES,
                    GATKSVVCFConstants.INSERTED_SEQUENCE_MAPPINGS,
                    GATKSVVCFConstants.TOTAL_MAPPINGS,
                    GATKSVVCFConstants.SPLIT_READ_SUPPORT,
                    GATKSVVCFConstants.READ_PAIR_SUPPORT);

    private static final class DiscoverVariantsFromContigAlignmentsSAMSparkIntegrationTestArgs {
        final String outputDir;
        final String cnvCallsLoc;

        DiscoverVariantsFromContigAlignmentsSAMSparkIntegrationTestArgs(final String outputDir, final String cnvCallsLoc){
            this.outputDir = outputDir;
            this.cnvCallsLoc = cnvCallsLoc;
        }

        String getCommandLine() {
            return  " -R " + SVIntegrationTestDataProvider.reference_2bit +
                    " -I " + SVIntegrationTestDataProvider.TEST_CONTIG_SAM +
                    " -O " + outputDir + "/DiscoverVariantsFromContigAlignmentsSAMSparkIntegrationTest" +
                    (cnvCallsLoc == null ? "" : " --cnv-calls " + cnvCallsLoc);
        }

    }

    @DataProvider(name = "discoverVariantsFromContigAlignmentsSparkIntegrationTest")
    public Object[][] createTestData() {
        List<Object[]> tests = new ArrayList<>();
        final File tempDirLeft = GATKBaseTest.createTempDir("forLeft");
        tempDirLeft.deleteOnExit();
        tests.add(new Object[]{
                new DiscoverVariantsFromContigAlignmentsSAMSparkIntegrationTestArgs(tempDirLeft.getAbsolutePath(), SVIntegrationTestDataProvider.EXTERNAL_CNV_CALLS)
        });
        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "discoverVariantsFromContigAlignmentsSparkIntegrationTest", groups = "sv")
    public void testDiscoverVariantsRunnableLocal(final DiscoverVariantsFromContigAlignmentsSAMSparkIntegrationTestArgs params) throws IOException {

        final List<String> args = Arrays.asList( new ArgumentsBuilder().add(params.getCommandLine()).getArgsArray() );
        runCommandLine(args);
        final String newVCF = args.get(args.indexOf("-O") + 1) + "_sample_inv_del_ins.vcf";
        StructuralVariationDiscoveryPipelineSparkIntegrationTest.svDiscoveryVCFEquivalenceTest(newVCF,
                SVIntegrationTestDataProvider.EXPECTED_SIMPLE_DEL_VCF, null, annotationsToIgnoreWhenComparingVariants, false);
    }

    @Test(dataProvider = "discoverVariantsFromContigAlignmentsSparkIntegrationTest", groups = "sv")
    public void testDiscoverVariantsRunnableMiniCluster(final DiscoverVariantsFromContigAlignmentsSAMSparkIntegrationTestArgs params) throws Exception {

        MiniClusterUtils.runOnIsolatedMiniCluster(cluster -> {

            final List<String> argsToBeModified = Arrays.asList( new ArgumentsBuilder().add(params.getCommandLine()).getArgsArray() );
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

            // outputs, prefix with hdfs address
            idx = argsToBeModified.indexOf("-O");
            path = new Path(workingDirectory, "test");
            final String vcfOnHDFS = path.toUri().toString() + "_sample_inv_del_ins.vcf";
            argsToBeModified.set(idx+1, path.toUri().toString());

            runCommandLine(argsToBeModified);
            StructuralVariationDiscoveryPipelineSparkIntegrationTest.svDiscoveryVCFEquivalenceTest(
                    vcfOnHDFS,
                    SVIntegrationTestDataProvider.EXPECTED_SIMPLE_DEL_VCF,
                    null,
                    annotationsToIgnoreWhenComparingVariants,
                    true);
        });
    }
}
