package org.broadinstitute.hellbender.tools.spark.sv.integration;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import org.apache.hadoop.fs.Path;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.test.MiniClusterUtils;
import org.broadinstitute.hellbender.utils.test.VariantContextTestUtils;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

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
    public void testSVDiscoverPipelineRunnableLocal(final StructuralVariationDiscoveryPipelineSparkIntegrationTest.StructuralVariationDiscoveryPipelineSparkIntegrationTestArgs params) throws Exception {

        final List<String> args = Arrays.asList( new ArgumentsBuilder().add(params.getCommandLineNoApiKey()).getArgsArray() );
        runCommandLine(args);

        svDiscoveryVCFEquivalenceTest(args.get(args.indexOf("-O")+1), SVIntegrationTestDataProvider.EXPECTED_SIMPLE_DEL_VCF, Arrays.asList("ALIGN_LENGTHS", "CTG_NAMES"), false);
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
            final String vcfOnHDFS = path.toUri().toString();
            argsToBeModified.set(idx+1, vcfOnHDFS);

            idx = argsToBeModified.indexOf("--contigSAMFile");
            path = new Path(workingDirectory, "assemblies.sam");
            argsToBeModified.set(idx+1, path.toUri().toString());

            idx = argsToBeModified.indexOf("--breakpointIntervals");
            path = new Path(workingDirectory, "intervals");
            argsToBeModified.set(idx+1, path.toUri().toString());

            idx = argsToBeModified.indexOf("--fastqDir");
            path = new Path(workingDirectory, "fastq");
            argsToBeModified.set(idx+1, path.toUri().toString());

            runCommandLine(argsToBeModified);
            svDiscoveryVCFEquivalenceTest(vcfOnHDFS, SVIntegrationTestDataProvider.EXPECTED_SIMPLE_DEL_VCF, Arrays.asList("ALIGN_LENGTHS", "CTG_NAMES"), true);
        });
    }

    public static void svDiscoveryVCFEquivalenceTest(final String generatedVCFPath, final String expectedVCFPath,
                                                     final List<String> attributesToIgnore, final boolean onHDFS) throws Exception{

        VCFFileReader fileReader;
        CloseableIterator<VariantContext> iterator;
        final List<VariantContext> actualVcs;
        if (onHDFS) {
            final File tempLocalVCF = BaseTest.createTempFile("variants", "vcf");
            tempLocalVCF.deleteOnExit();
            BucketUtils.copyFile(generatedVCFPath, tempLocalVCF.getAbsolutePath());
            fileReader = new VCFFileReader(tempLocalVCF, false);
        } else {
            fileReader = new VCFFileReader(new File(generatedVCFPath), false);
        }
        iterator = fileReader.iterator();
        actualVcs = Utils.stream(iterator).collect(Collectors.toList());
        CloserUtil.close(iterator);
        CloserUtil.close(fileReader);

        fileReader = new VCFFileReader(new File(expectedVCFPath), false);
        iterator = fileReader.iterator();
        final List<VariantContext> expectedVcs = Utils.stream(iterator).collect(Collectors.toList());
        CloserUtil.close(iterator);
        CloserUtil.close(fileReader);

        BaseTest.assertCondition(actualVcs, expectedVcs, (a, e) -> VariantContextTestUtils.assertVariantContextsAreEqual(a, e, attributesToIgnore));
    }
}
