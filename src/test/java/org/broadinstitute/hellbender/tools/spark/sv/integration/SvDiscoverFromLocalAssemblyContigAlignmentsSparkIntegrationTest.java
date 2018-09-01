package org.broadinstitute.hellbender.tools.spark.sv.integration;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SvDiscoverFromLocalAssemblyContigAlignmentsSpark;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class SvDiscoverFromLocalAssemblyContigAlignmentsSparkIntegrationTest extends CommandLineProgramTest {

    private static final String REFERENCE = largeFileTestDir + "/Homo_sapiens_assembly38.20.21.2bit";
    private static final String NON_CANONICAL_CHROMOSOMES = toolsTestDir + "spark/sv/integration/inputs/Homo_sapiens_assembly38.kill.alts";

    private static final String THIS_TEST_DIR = largeFileTestDir + "/sv/";
    private static final String TEST_SAMPLE_NAME = "SvDiscoverFromLocalAssemblyContigAlignmentsSparkIntegrationTestSample";

    // TODO: 8/15/18 this split between cases of BND variant and "symbolic" variants is temporary because of restraint placed by access to only chr20 and chr21 (BND variant are lacking cases); in the long run we'd prefer to merge them
    private static final String SYMBOLIC_ASSEMBLY_BAM = THIS_TEST_DIR + "SvDiscoverFromLocalAssemblyContigAlignmentsSparkIntegrationTest.bam";
    private static final String SYMBOLIC_EXPECTED_SIMPLE_CHIMERA_VCF = THIS_TEST_DIR + "SvDiscoverFromLocalAssemblyContigAlignmentsSparkIntegrationTest_NonComplex.vcf";
    private static final String SYMBOLIC_EXPECTED_COMPLEX_CHIMERA_VCF = THIS_TEST_DIR + "SvDiscoverFromLocalAssemblyContigAlignmentsSparkIntegrationTest_Complex.vcf";
    private static final String SYMBOLIC_EXPECTED_CPX_RETERINTERPRETATION_1_SEG_VCF = THIS_TEST_DIR + "SvDiscoverFromLocalAssemblyContigAlignmentsSparkIntegrationTest_cpx_reinterpreted_simple_1_seg.vcf";
    private static final String SYMBOLIC_EXPECTED_CPX_RETERINTERPRETATION_MULTI_SEG_VCF = THIS_TEST_DIR + "SvDiscoverFromLocalAssemblyContigAlignmentsSparkIntegrationTest_cpx_reinterpreted_simple_multi_seg.vcf";
    private static final String SYMBOLIC_EXPECTED_MERGED_SIMPLE_VCF = THIS_TEST_DIR + "SvDiscoverFromLocalAssemblyContigAlignmentsSparkIntegrationTest_merged_simple.vcf";

    private static final String BND_ASSEMBLY_BAM = THIS_TEST_DIR + "forBND_SvDiscoverFromLocalAssemblyContigAlignmentsSparkIntegrationTest.bam";
    private static final String BND_EXPECTED_SIMPLE_CHIMERA_VCF = THIS_TEST_DIR + "forBND_SvDiscoverFromLocalAssemblyContigAlignmentsSparkIntegrationTest_NonComplex.vcf";
    private static final String BND_EXPECTED_MERGED_SIMPLE_VCF = THIS_TEST_DIR + "forBND_SvDiscoverFromLocalAssemblyContigAlignmentsSparkIntegrationTest_merged_simple.vcf";

    private static final class TestArgs {

        final String outputDir;
        final String twoBitReference;
        final String altContigsFile;
        final String assemblyBam;
        final String cnvCallsLoc;

        TestArgs(final String workingDir, final String twoBitReference, final String altContigsFile, final String assemblyBam, final String cnvCallsLoc) {
            this.outputDir = workingDir.endsWith("/") ? workingDir : workingDir + "/";
            this.twoBitReference = twoBitReference;
            this.altContigsFile = altContigsFile;
            this.assemblyBam = assemblyBam;
            this.cnvCallsLoc = cnvCallsLoc;
        }

        String getCommandLine() {
            return  " --debug-mode " +
                    " -R " + twoBitReference +
                    " -I " + assemblyBam +
                    " -O " + outputDir +
                    (altContigsFile == null ? "" : " -alt-tigs " + altContigsFile) +
                    (cnvCallsLoc == null ? "" : " --cnv-calls " + cnvCallsLoc);
        }

        String getOutputDir() {
            return outputDir;
        }
    }

    @DataProvider(name = "svDiscoverFromLocalAssemblyContigAlignmentsSparkIntegrationTest")
    public Object[][] createTestData() {
        final File tempDir = GATKBaseTest.createTempDir("symb");
        tempDir.deleteOnExit();

        final List<Object[]> tests = new ArrayList<>(2);
        final TestArgs testArgs = new TestArgs(tempDir.getAbsolutePath(), REFERENCE, NON_CANONICAL_CHROMOSOMES,
                SYMBOLIC_ASSEMBLY_BAM, null);
        tests.add(new Object[]{testArgs,
                SYMBOLIC_EXPECTED_SIMPLE_CHIMERA_VCF,
                SYMBOLIC_EXPECTED_COMPLEX_CHIMERA_VCF,
                SYMBOLIC_EXPECTED_CPX_RETERINTERPRETATION_1_SEG_VCF,
                SYMBOLIC_EXPECTED_CPX_RETERINTERPRETATION_MULTI_SEG_VCF,
                SYMBOLIC_EXPECTED_MERGED_SIMPLE_VCF
        });

        final File tempDir1 = GATKBaseTest.createTempDir("bnd");
        tempDir1.deleteOnExit();
        final TestArgs testArgs1 = new TestArgs(tempDir1.getAbsolutePath(), REFERENCE, NON_CANONICAL_CHROMOSOMES,
                BND_ASSEMBLY_BAM, null);
        tests.add(new Object[]{testArgs1,
                BND_EXPECTED_SIMPLE_CHIMERA_VCF,
                null,
                null,
                null,
                BND_EXPECTED_MERGED_SIMPLE_VCF
        });

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "svDiscoverFromLocalAssemblyContigAlignmentsSparkIntegrationTest", groups = "sv")
    public void testDiscoverVariantsRunnableLocal(final TestArgs params, final String expectedSimpleChimeraFilePath,
                                                  final String expectedComplexFilePath, final String expectedReInterpretationOneSegFilePath, final String expectedReInterPretationMultiSegFilePath,
                                                  final String expectedMergedSimpleFilePath) throws IOException {

        final List<String> args = Arrays.asList( new ArgumentsBuilder().add(params.getCommandLine()).getArgsArray() );
        runCommandLine(args);

        final String outputDir = params.getOutputDir();

        final String simpleChimeraVCF = outputDir + TEST_SAMPLE_NAME + "_" + SvDiscoverFromLocalAssemblyContigAlignmentsSpark.SIMPLE_CHIMERA_VCF_FILE_NAME;
        StructuralVariationDiscoveryPipelineSparkIntegrationTest
                .svDiscoveryVCFEquivalenceTest(
                        simpleChimeraVCF,
                        expectedSimpleChimeraFilePath,
                        null,
                        Collections.emptyList(),
                        false);

        final String complexChimeraVCF = outputDir + TEST_SAMPLE_NAME + "_" + SvDiscoverFromLocalAssemblyContigAlignmentsSpark.COMPLEX_CHIMERA_VCF_FILE_NAME;
        StructuralVariationDiscoveryPipelineSparkIntegrationTest
                .svDiscoveryVCFEquivalenceTest(
                        complexChimeraVCF,
                        expectedComplexFilePath,
                        null,
                        Collections.emptyList(),
                        false);

        final String cpxReInterpretOneSegVCF = outputDir + TEST_SAMPLE_NAME + "_" + SvDiscoverFromLocalAssemblyContigAlignmentsSpark.REINTERPRETED_1_SEG_CALL_VCF_FILE_NAME;
        StructuralVariationDiscoveryPipelineSparkIntegrationTest
                .svDiscoveryVCFEquivalenceTest(
                        cpxReInterpretOneSegVCF,
                        expectedReInterpretationOneSegFilePath,
                        null,
                        Collections.emptyList(),
                        false);

        final String cpxReInterpretMultiSegVCF = outputDir + TEST_SAMPLE_NAME + "_" + SvDiscoverFromLocalAssemblyContigAlignmentsSpark.REINTERPRETED_MULTI_SEG_CALL_VCF_FILE_NAME;
        StructuralVariationDiscoveryPipelineSparkIntegrationTest
                .svDiscoveryVCFEquivalenceTest(
                        cpxReInterpretMultiSegVCF,
                        expectedReInterPretationMultiSegFilePath,
                        null,
                        Collections.emptyList(),
                        false);

        final String mergedSimpleVCF = outputDir + TEST_SAMPLE_NAME + "_" + SvDiscoverFromLocalAssemblyContigAlignmentsSpark.MERGED_VCF_FILE_NAME;
        StructuralVariationDiscoveryPipelineSparkIntegrationTest
                .svDiscoveryVCFEquivalenceTest(
                        mergedSimpleVCF,
                        expectedMergedSimpleFilePath,
                        null,
                        Collections.emptyList(),
                        false);
    }
}
