package org.broadinstitute.hellbender.tools.spark.sv.integration;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import org.apache.hadoop.fs.Path;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.MiniClusterUtils;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

public class CpxVariantReInterpreterSparkIntegrationTest extends CommandLineProgramTest {

    private static final String THIS_TEST_FOLDER = getTestDataDir() + "/spark/sv/integration/inputs/";

    private static final String complexVCF = THIS_TEST_FOLDER + "CpxVariantReInterpreterSparkIntegrationTest_complex.vcf";
    private static final String assemblyAlignmentsAccompanyingComplexVCF = largeFileTestDir + "/sv/CpxVariantReInterpreterSparkIntegrationTest_complex_assembly.bam";
    private static final File fastaReference = new File(b38_reference_20_21);
    private static final String nonCanonicalChromosomeNamesFile = THIS_TEST_FOLDER + "Homo_sapiens_assembly38.kill.alts";
    private static final String outPrefix = "cpx_reinterpreted_simple";
    private static final String expectedOneSegmentVCF = getTestDataDir() + "/spark/sv/integration/outputs/cpx_reinterpreted_simple_1_seg.vcf";
    private static final String expectedMultiSegmentVCF = getTestDataDir() + "/spark/sv/integration/outputs/cpx_reinterpreted_simple_multi_seg.vcf";
    private static final List<String> annotationsToIgnoreWhenComparingVariants = Collections.emptyList();

    private static final class CpxVariantReInterpreterSparkIntegrationTestArgs {
        final String outputDir;

        CpxVariantReInterpreterSparkIntegrationTestArgs(final String outputDir) {
            this.outputDir = outputDir;
        }

        String getCommandLine() {
            return  " -R " + fastaReference +
                    " -I " + assemblyAlignmentsAccompanyingComplexVCF +
                    " --non-canonical-contig-names-file " + nonCanonicalChromosomeNamesFile +
                    " --cpx-vcf " + complexVCF +
                    " --prefix-out-vcf " + outputDir + "/" + outPrefix;
        }
    }

    @DataProvider(name = "forCpxVariantReInterpreterSparkIntegrationTest")
    private Object[][] createData() {
        List<Object[]> data = new ArrayList<>();
        final File cpxVariantReInterpreterSparkIntegrationTest = createTempDir("CpxVariantReInterpreterSparkIntegrationTest");
        data.add(new Object[]{new CpxVariantReInterpreterSparkIntegrationTestArgs(cpxVariantReInterpreterSparkIntegrationTest.getAbsolutePath())});
        return data.toArray(new Object[data.size()][]);
    }

    @Test(groups = "sv", dataProvider = "forCpxVariantReInterpreterSparkIntegrationTest")
    public void testRunLocal(final CpxVariantReInterpreterSparkIntegrationTestArgs params) throws Exception {
        final List<String> args = Arrays.asList( new ArgumentsBuilder().add(params.getCommandLine()).getArgsArray() );
        runCommandLine(args);

        final String actualVCFForOneSegmentCalls = params.outputDir + "/" + outPrefix + "_1_seg.vcf";
        vcfEquivalenceTest(actualVCFForOneSegmentCalls, expectedOneSegmentVCF, annotationsToIgnoreWhenComparingVariants, false);
        final String actualVCFForMultiSegmentCalls = params.outputDir + "/" + outPrefix + "_multi_seg.vcf";
        vcfEquivalenceTest(actualVCFForMultiSegmentCalls, expectedMultiSegmentVCF, annotationsToIgnoreWhenComparingVariants, false);
    }

    @Test(groups = "sv", dataProvider = "forCpxVariantReInterpreterSparkIntegrationTest")
    public void testRunHDFS(final CpxVariantReInterpreterSparkIntegrationTestArgs params) throws Exception {
        MiniClusterUtils.runOnIsolatedMiniCluster(cluster -> {

            final List<String> argsToBeModified = Arrays.asList( new ArgumentsBuilder().add(params.getCommandLine()).getArgsArray() );
            final Path workingDirectory = MiniClusterUtils.getWorkingDir(cluster);

            int idx = 0;

            // copy inputs
            idx = argsToBeModified.indexOf("-I");
            Path path = new Path(workingDirectory, "hdfs.bam");
            File file = new File(argsToBeModified.get(idx+1));
            cluster.getFileSystem().copyFromLocalFile(new Path(file.toURI()), path);
            argsToBeModified.set(idx+1, path.toUri().toString());
            path = new Path(workingDirectory, "hdfs.bam.bai"); // .bai
            file = new File(file.getAbsolutePath() + ".bai");
            cluster.getFileSystem().copyFromLocalFile(new Path(file.toURI()), path);

            idx = argsToBeModified.indexOf("-R");
            path = new Path(workingDirectory, "reference.fasta");
            file = new File(argsToBeModified.get(idx+1));
            cluster.getFileSystem().copyFromLocalFile(new Path(file.toURI()), path);
            argsToBeModified.set(idx+1, path.toUri().toString());
            path = new Path(workingDirectory, "reference.fasta.fai"); // .fasta.fai for fasta
            file = new File(file.getAbsolutePath() + ".fai");
            cluster.getFileSystem().copyFromLocalFile(new Path(file.toURI()), path);
            path = new Path(workingDirectory, "reference.dict"); // .dict for fasta
            file = new File(file.getAbsolutePath().replace(".fasta.fai", ".dict"));
            cluster.getFileSystem().copyFromLocalFile(new Path(file.toURI()), path);

            idx = argsToBeModified.indexOf("--non-canonical-contig-names-file");
            path = new Path(workingDirectory, "Homo_sapiens_assembly38.kill.alts");
            file = new File(argsToBeModified.get(idx+1));
            cluster.getFileSystem().copyFromLocalFile(new Path(file.toURI()), path);
            argsToBeModified.set(idx+1, path.toUri().toString());

            idx = argsToBeModified.indexOf("--cpx-vcf");
            path = new Path(workingDirectory, "CpxVariantReInterpreterSparkIntegrationTest_complex.vcf");
            file = new File(argsToBeModified.get(idx+1));
            cluster.getFileSystem().copyFromLocalFile(new Path(file.toURI()), path);
            argsToBeModified.set(idx+1, path.toUri().toString());

            // outputs, prefix with hdfs address
            idx = argsToBeModified.indexOf("--prefix-out-vcf");
            path = new Path(workingDirectory, "test");
            argsToBeModified.set(idx+1, path.toUri().toString());

            runCommandLine(argsToBeModified);

            final String actualVCFForOneSegmentCallsOnHDFS = path.toUri().toString() + "_1_seg.vcf";
            final String actualVCFForMultiSegmentCallsOnHDFS = path.toUri().toString() + "_multi_seg.vcf";

            vcfEquivalenceTest(actualVCFForOneSegmentCallsOnHDFS, expectedOneSegmentVCF, annotationsToIgnoreWhenComparingVariants, true);
            vcfEquivalenceTest(actualVCFForMultiSegmentCallsOnHDFS, expectedMultiSegmentVCF, annotationsToIgnoreWhenComparingVariants, true);
        });
    }

    private static void vcfEquivalenceTest(final String generatedVCFPath, final String expectedVCFPath,
                                           final List<String> attributesToIgnore, final boolean onHDFS) throws Exception {

        List<VariantContext> expectedVcs;
        try (final VCFFileReader fileReader = new VCFFileReader(new File(expectedVCFPath), false) ) {
            try (final CloseableIterator<VariantContext> iterator = fileReader.iterator()) {
                expectedVcs = Utils.stream(iterator).collect(Collectors.toList());
            }
        }

        final List<VariantContext> actualVcs = StructuralVariationDiscoveryPipelineSparkIntegrationTest
                .extractActualVCs(generatedVCFPath, onHDFS);

        GATKBaseTest.assertCondition(actualVcs, expectedVcs,
                (a, e) -> VariantContextTestUtils.assertVariantContextsAreEqual(a, e, attributesToIgnore));
    }
}
