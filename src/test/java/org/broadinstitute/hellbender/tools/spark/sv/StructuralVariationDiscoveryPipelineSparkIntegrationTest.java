package org.broadinstitute.hellbender.tools.spark.sv;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.test.IntegrationTestSpec;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
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
    public void testSVDiscoverPipelineRunnableMiniCluster(final StructuralVariationDiscoveryPipelineSparkIntegrationTest.StructuralVariationDiscoveryPipelineSparkIntegrationTestArgs params) throws IOException {

    }
}
