package org.broadinstitute.hellbender.tools.spark.sv.integration;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import org.apache.hadoop.fs.Path;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
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
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

import static org.broadinstitute.hellbender.tools.spark.sv.integration.DiscoverVariantsFromContigAlignmentsSAMSparkIntegrationTest.annotationsToIgnoreWhenComparingVariants;

public class StructuralVariationDiscoveryPipelineSparkIntegrationTest extends CommandLineProgramTest {

    private static final class StructuralVariationDiscoveryPipelineSparkIntegrationTestArgs {
        private static final String EXPERIMENTAL_INTERPRETATION_OUTPUT_DIR_NAME = "experimentalVariantInterpretations";

        final String bamLoc;
        final String kmerIgnoreListLoc;
        final String alignerRefIndexImgLoc;
        final String outputDir;
        final String cnvCallsLoc;


        StructuralVariationDiscoveryPipelineSparkIntegrationTestArgs(final String bamLoc,
                                                                     final String kmerIgnoreListLoc,
                                                                     final String alignerRefIndexImgLoc,
                                                                     final String cnvCallsLoc,
                                                                     final String outputDir) {
            this.bamLoc = bamLoc;
            this.kmerIgnoreListLoc = kmerIgnoreListLoc;
            this.alignerRefIndexImgLoc = alignerRefIndexImgLoc;
            this.outputDir = outputDir;
            this.cnvCallsLoc = cnvCallsLoc;
        }

        String getCommandLine() {
            return  " -R " + SVIntegrationTestDataProvider.reference_2bit +
                    " -I " + bamLoc +
                    " -O " + outputDir        + "/variants.vcf" +
            " --aligner-index-image " + alignerRefIndexImgLoc +
                    " --kmers-to-ignore " + kmerIgnoreListLoc +
                    " --contig-sam-file "       + outputDir + "/assemblies.sam" +
                    " --breakpoint-intervals " + outputDir + "/intervals" +
                    " --fastq-dir "            + outputDir + "/fastq" +
                    (cnvCallsLoc == null ? "" : " --cnv-calls " + cnvCallsLoc) +
                    " --exp-variants-out-dir " + getExperimentalInterpretationOutputDirName(outputDir + "/variants.vcf");
        }

        private static String getExperimentalInterpretationOutputDirName(final String vcfPath) {
            return Paths.get(vcfPath).getParent().resolve(EXPERIMENTAL_INTERPRETATION_OUTPUT_DIR_NAME)
                    .toAbsolutePath().toString();
        }

        @Override
        public String toString() {
            return "StructuralVariationDiscoveryPipelineSparkIntegrationTestArgs{" +
                    "bam-loc='" + bamLoc + '\'' +
                    ", kmer-ignore-list-loc='" + kmerIgnoreListLoc + '\'' +
                    ", aligner-fef-index-img-loc='" + alignerRefIndexImgLoc + '\'' +
                    ", cnv-calls-loc='" + cnvCallsLoc + '\'' +
                    ", output-dir='" + outputDir + '\'' +
                    '}';
        }
    }

    @DataProvider(name = "svDiscoverPipelineSparkIntegrationTest")
    public Object[][] createTestData() throws IOException {
        List<Object[]> tests = new ArrayList<>();

        final File tempDirNew = BaseTest.createTempDir("new");
        tempDirNew.deleteOnExit();
        Files.createDirectories(Paths.get(tempDirNew.getAbsolutePath()+"/fastq"));
        tests.add(new Object[]{
                new StructuralVariationDiscoveryPipelineSparkIntegrationTest.StructuralVariationDiscoveryPipelineSparkIntegrationTestArgs(
                        SVIntegrationTestDataProvider.TEST_BAM,
                        SVIntegrationTestDataProvider.KMER_KILL_LIST,
                        SVIntegrationTestDataProvider.ALIGNER_INDEX_IMG,
                        SVIntegrationTestDataProvider.EXTERNAL_CNV_CALLS,
                        tempDirNew.getAbsolutePath()
                )
        });

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "svDiscoverPipelineSparkIntegrationTest", groups = "sv")
    public void testSVDiscoverPipelineRunnableLocal(final StructuralVariationDiscoveryPipelineSparkIntegrationTest.StructuralVariationDiscoveryPipelineSparkIntegrationTestArgs params) throws Exception {

        final List<String> args = Arrays.asList( new ArgumentsBuilder().add(params.getCommandLine()).getArgsArray() );
        runCommandLine(args);

        svDiscoveryVCFEquivalenceTest(args.get(args.indexOf("-O")+1), SVIntegrationTestDataProvider.EXPECTED_SIMPLE_DEL_VCF,
                args.get(args.indexOf("--exp-variants-out-dir")+1),
                annotationsToIgnoreWhenComparingVariants, false);
    }

    @Test(dataProvider = "svDiscoverPipelineSparkIntegrationTest", groups = "sv")
    public void testSVDiscoverPipelineRunnableMiniCluster(final StructuralVariationDiscoveryPipelineSparkIntegrationTest.StructuralVariationDiscoveryPipelineSparkIntegrationTestArgs params) throws Exception {

        MiniClusterUtils.runOnIsolatedMiniCluster(cluster -> {

            final List<String> argsToBeModified = Arrays.asList( new ArgumentsBuilder().add(params.getCommandLine()).getArgsArray() );
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

            idx = argsToBeModified.indexOf("--kmers-to-ignore");
            path = new Path(workingDirectory, "dummy.kill.kmers");
            file = new File(argsToBeModified.get(idx+1));
            cluster.getFileSystem().copyFromLocalFile(new Path(file.toURI()), path);
            argsToBeModified.set(idx+1, path.toUri().toString());

            idx = argsToBeModified.indexOf("--cnv-calls");
            path = new Path(workingDirectory, "cnvVariants");
            file = new File(argsToBeModified.get(idx+1));
            cluster.getFileSystem().copyFromLocalFile(new Path(file.toURI()), path);
            argsToBeModified.set(idx+1, path.toUri().toString());

            // outputs, prefix with hdfs address
            idx = argsToBeModified.indexOf("-O");
            path = new Path(workingDirectory, "variants.vcf");
            final String vcfOnHDFS = path.toUri().toString();
            argsToBeModified.set(idx+1, vcfOnHDFS);

            idx = argsToBeModified.indexOf("--contig-sam-file");
            path = new Path(workingDirectory, "assemblies.sam");
            argsToBeModified.set(idx+1, path.toUri().toString());

            idx = argsToBeModified.indexOf("--breakpoint-intervals");
            path = new Path(workingDirectory, "intervals");
            argsToBeModified.set(idx+1, path.toUri().toString());

            idx = argsToBeModified.indexOf("--fastq-dir");
            path = new Path(workingDirectory, "fastq");
            argsToBeModified.set(idx+1, path.toUri().toString());

            idx = argsToBeModified.indexOf("--exp-variants-out-dir");
            path = new Path(workingDirectory, "expVariantsOutDir");
            final String expOutDirOnHDFS = path.toUri().toString();
            argsToBeModified.set(idx+1, expOutDirOnHDFS);

            runCommandLine(argsToBeModified);
            svDiscoveryVCFEquivalenceTest(vcfOnHDFS, SVIntegrationTestDataProvider.EXPECTED_SIMPLE_DEL_VCF,
                    expOutDirOnHDFS, annotationsToIgnoreWhenComparingVariants, true);
        });
    }

    static void svDiscoveryVCFEquivalenceTest(final String generatedVCFPath, final String expectedVCFPath,
                                              final String experimentalOutputPath,
                                              final List<String> attributesToIgnore, final boolean onHDFS) throws Exception {

        final VCFFileReader fileReader = new VCFFileReader(new File(expectedVCFPath), false);
        final CloseableIterator<VariantContext> iterator = fileReader.iterator();
        final List<VariantContext> expectedVcs = Utils.stream(iterator).collect(Collectors.toList());
        CloserUtil.close(iterator);
        CloserUtil.close(fileReader);

        List<VariantContext> actualVcs = extractActualVCs(generatedVCFPath, onHDFS);

        GATKBaseTest.assertCondition(actualVcs, expectedVcs,
                (a, e) -> VariantContextTestUtils.assertVariantContextsAreEqual(a, e, attributesToIgnore));

        if ( experimentalOutputPath != null ) {
            final java.nio.file.Path path = IOUtils.getPath(experimentalOutputPath).resolve("nonComplex.vcf");
            final String experimentalInsDelVcf = onHDFS ? path.toUri().toString() : path.toString();
            actualVcs = extractActualVCs(experimentalInsDelVcf, onHDFS);

            // TODO: 1/28/18 see ticket #4228
            final List<String> moreAttributesToIgnoreForNow = new ArrayList<>(attributesToIgnore);
            moreAttributesToIgnoreForNow.addAll(Collections.singletonList("EXTERNAL_CNV_CALLS"));
            GATKBaseTest.assertCondition(actualVcs, expectedVcs,
                    (a, e) -> VariantContextTestUtils.assertVariantContextsAreEqual(a, e, moreAttributesToIgnoreForNow));
        }
    }

    private static List<VariantContext> extractActualVCs(final String generatedVCFPath, final boolean onHDFS)
            throws IOException {

        final VCFFileReader fileReader;
        if (onHDFS) {
            final File tempLocalVCF = GATKBaseTest.createTempFile("variants", "vcf");
            tempLocalVCF.deleteOnExit();
            BucketUtils.copyFile(generatedVCFPath, tempLocalVCF.getAbsolutePath());
            fileReader = new VCFFileReader(tempLocalVCF, false);
        } else {
            fileReader = new VCFFileReader(new File(generatedVCFPath), false);
        }
        final CloseableIterator<VariantContext> iterator = fileReader.iterator();
        final List<VariantContext> actualVcs = Utils.stream(iterator).collect(Collectors.toList());
        CloserUtil.close(iterator);
        CloserUtil.close(fileReader);

        return actualVcs;
    }
}
