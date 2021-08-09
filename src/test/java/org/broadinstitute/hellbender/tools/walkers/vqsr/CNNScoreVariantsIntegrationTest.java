package org.broadinstitute.hellbender.tools.walkers.vqsr;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
import org.broadinstitute.hellbender.utils.Utils;

import org.broadinstitute.hellbender.utils.python.PythonScriptExecutorException;
import org.testng.Assert;
import org.testng.SkipException;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

import java.util.Collections;
import java.util.Iterator;
import java.util.stream.Collectors;


/**
 * Integration tests for {@link CNNScoreVariants}.
 * Created by sam on 1/8/18.
 */
public class CNNScoreVariantsIntegrationTest extends CommandLineProgramTest {
    private static final String architecture1D = largeFileTestDir + "VQSR/cnn_ref_model/1d_cnn_mix_train_full_bn.json";
    private static final String weights1D = largeFileTestDir + "VQSR/cnn_ref_model/1d_cnn_mix_train_full_bn.hd5";
    private static final String architecture2D = largeFileTestDir + "VQSR/cnn_read_model/small_2d.json";
    private static final String weights2D = largeFileTestDir + "VQSR/cnn_read_model/small_2d.hd5";
    private static final String inputVCF = largeFileTestDir + "VQSR/recalibrated_chr20_start.vcf";
    private static final String bigInputVCF = largeFileTestDir + "VQSR/g94982_20_1m_10m_python_2dcnn.vcf.gz";
    private static final String inputBAM = largeFileTestDir + "VQSR/g94982_contig_20_start_bamout.bam";
    private static final String inputIntervals = largeFileTestDir + "VQSR/contig20_conf_1m_10m.interval_list";
    private static final double EPSILON = 0.01;

    @Test(expectedExceptions = RuntimeException.class)
    public void testRequirePythonEnvironment() throws IOException {
        // This test is deliberately left out of the "python" test group in order to ensure that
        // it only executes when the Python environment has *NOT* been properly established. Also,
        // skip this test if we're running on the Docker because the Python environment is always
        // activated there.
        if (isGATKDockerContainer()) {
            throw new SkipException("Python environment validation test must be skipped when running on the Docker");
        }

        // Re-running the "testAllDefaultArgs" test should throw when run outside of the GATK Python environment
        testAllDefaultArgs();
    }
    /**
     * Run the tool on a small test VCF.
     */
    @Test(groups = {"python"})
    public void testAllDefaultArgs() {
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        final File tempVcf = createTempFile("tester", ".vcf");
        final File expectedVcf = new File(largeFileTestDir + "VQSR/expected/cnn_1d_chr20_subset_expected.vcf");
        argsBuilder.add(StandardArgumentDefinitions.VARIANT_LONG_NAME, inputVCF)
                .add(StandardArgumentDefinitions.OUTPUT_LONG_NAME, tempVcf.getPath())
                .add(StandardArgumentDefinitions.REFERENCE_LONG_NAME, b37_reference_20_21)
                .add(StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false");

        runCommandLine(argsBuilder);
        assertInfoFieldsAreClose(tempVcf, expectedVcf, GATKVCFConstants.CNN_1D_KEY);
    }

    @Test(groups = {"python"}, expectedExceptions = PythonScriptExecutorException.class)
    public void testExceptionDuringAsyncBatch() {
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        final File tempVcf = createTempFile("tester", ".vcf");
        // the last variant in this vcf has a  value of "." for the float attributes in the default CNN
        // annotation set MQ, MQRankSum, ReadPosRankSum, SOR, VQSLOD, and QD
        //TODO: move this into the large resources dir
        final File malformedVCF = new File("src/test/resources/cnn_1d_chr20_subset_expected.badAnnotations.vcf");
        argsBuilder.add(StandardArgumentDefinitions.VARIANT_LONG_NAME, malformedVCF)
                .add(StandardArgumentDefinitions.OUTPUT_LONG_NAME, tempVcf.getPath())
                .add(StandardArgumentDefinitions.REFERENCE_LONG_NAME, b37_reference_20_21)
                .add(StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false");

        runCommandLine(argsBuilder);
    }

    @Test(groups = {"python"})
    public void testInferenceArchitecture() {
        final boolean newExpectations = false;
        final String expectedVCFName = largeFileTestDir + "VQSR/expected/cnn_1d_chr20_subset_expected.vcf";
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        argsBuilder.add(StandardArgumentDefinitions.VARIANT_LONG_NAME, inputVCF)
                .add(StandardArgumentDefinitions.REFERENCE_LONG_NAME, b37_reference_20_21)
                .add("architecture", architecture1D)
                .add(StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false");

        if (newExpectations) {
            argsBuilder.add(StandardArgumentDefinitions.OUTPUT_LONG_NAME, expectedVCFName);
            runCommandLine(argsBuilder);
        } else {
            final File tempVcf = createTempFile("tester", ".vcf");
            final File expectedVcf = new File(expectedVCFName);
            argsBuilder.add(StandardArgumentDefinitions.OUTPUT_LONG_NAME, tempVcf.getPath());
            runCommandLine(argsBuilder);
            assertInfoFieldsAreClose(tempVcf, expectedVcf, GATKVCFConstants.CNN_1D_KEY);
        }
    }

    @Test(groups = {"python"})
    public void testInferenceWeights() {
        final File tempVcf = createTempFile("tester", ".vcf");
        final File expectedVcf = new File(largeFileTestDir + "VQSR/expected/cnn_1d_chr20_subset_expected.vcf");
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        argsBuilder.add(StandardArgumentDefinitions.VARIANT_LONG_NAME, inputVCF)
                .add(StandardArgumentDefinitions.REFERENCE_LONG_NAME, b37_reference_20_21)
                .add("weights", weights1D)
                .add(StandardArgumentDefinitions.OUTPUT_LONG_NAME, tempVcf.getPath())
                .add(StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false");

        runCommandLine(argsBuilder);
        assertInfoFieldsAreClose(tempVcf, expectedVcf, GATKVCFConstants.CNN_1D_KEY);
    }

    @Test(groups = {"python"})
    public void testInferenceArchitectureAndWeights() {
        final File tempVcf = createTempFile("tester", ".vcf");
        final File expectedVcf = new File(largeFileTestDir + "VQSR/expected/cnn_1d_chr20_subset_expected.vcf");
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        argsBuilder.add(StandardArgumentDefinitions.VARIANT_LONG_NAME, inputVCF)
                .add(StandardArgumentDefinitions.REFERENCE_LONG_NAME, b37_reference_20_21)
                .add("weights", weights1D)
                .add("architecture", architecture1D)
                .add(StandardArgumentDefinitions.OUTPUT_LONG_NAME, tempVcf.getPath())
                .add(StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false");

        runCommandLine(argsBuilder);
        assertInfoFieldsAreClose(tempVcf, expectedVcf, GATKVCFConstants.CNN_1D_KEY);
    }

    @Test(groups = {"python"})
    public void testInferenceWithIntervals() {
        final boolean newExpectations = false;
        final String expectedVCFName = largeFileTestDir + "VQSR/expected/cnn_1d_contig20_1m_10m_expected.vcf";
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        argsBuilder.add(StandardArgumentDefinitions.VARIANT_LONG_NAME, bigInputVCF)
                .add(StandardArgumentDefinitions.REFERENCE_LONG_NAME, b37_reference_20_21)
                .add(StandardArgumentDefinitions.INTERVALS_LONG_NAME, inputIntervals)
                .add(StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false");

        if (newExpectations) {
            argsBuilder.add(StandardArgumentDefinitions.OUTPUT_LONG_NAME, expectedVCFName);
            runCommandLine(argsBuilder);
        } else {
            final File expectedVcf = new File(expectedVCFName);
            final File tempVcf = createTempFile("tester", ".vcf");
            argsBuilder.add(StandardArgumentDefinitions.OUTPUT_LONG_NAME, tempVcf.getPath());
            runCommandLine(argsBuilder);
            assertInfoFieldsAreClose(tempVcf, expectedVcf, GATKVCFConstants.CNN_1D_KEY);
        }
    }

    @Test(groups = {"python"})
    public void testSmallBatchInference() {
        final File tempVcf = createTempFile("tester", ".vcf");
        final File expectedVcf = new File(largeFileTestDir + "VQSR/expected/cnn_1d_chr20_subset_expected.vcf");
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        argsBuilder.add(StandardArgumentDefinitions.VARIANT_LONG_NAME, inputVCF)
                .add(StandardArgumentDefinitions.OUTPUT_LONG_NAME, tempVcf.getPath())
                .add(StandardArgumentDefinitions.REFERENCE_LONG_NAME, b37_reference_20_21)
                .add("inference-batch-size", "8")
                .add("transfer-batch-size", "16")
                .add(StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false");
        runCommandLine(argsBuilder);
        assertInfoFieldsAreClose(tempVcf, expectedVcf, GATKVCFConstants.CNN_1D_KEY);
    }

    @Test(groups = {"python"})
    public void testOnContigEdge() {
        final String edgeVcf = toolsTestDir + "walkers/VQSR/variantNearContigEdge.vcf";
        final File tempVcf = createTempFile("tester", ".vcf");
        final File expectedVcf = new File(largeFileTestDir + "VQSR/expected/chrM.vcf");
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        argsBuilder.add(StandardArgumentDefinitions.VARIANT_LONG_NAME, edgeVcf)
                .add(StandardArgumentDefinitions.OUTPUT_LONG_NAME, tempVcf.getPath())
                .add(StandardArgumentDefinitions.REFERENCE_LONG_NAME, hg19MiniReference)
                .add(StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false");

        runCommandLine(argsBuilder);
        assertInfoFieldsAreClose(tempVcf, expectedVcf, GATKVCFConstants.CNN_1D_KEY);
    }

    /**
     * Run the 2D Model on a small test VCF with the resource loaded weights and architecture.
     */
    @Test(groups = {"python"})
    public void testInference2dResourceModel() {
        // We reset the random number generator at the beginning of each test so that the random down-sampling of reads
        // by the reservoir down-sampler does not cause slightly different scores.
        Utils.resetRandomGenerator();
        TensorType tt = TensorType.read_tensor;
        final File tempVcf = createTempFile("tester", ".vcf");
        final File expectedVcf = new File(largeFileTestDir + "VQSR/expected/cnn_2d_chr20_subset_expected.vcf");
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        argsBuilder.add(StandardArgumentDefinitions.VARIANT_LONG_NAME, inputVCF)
                .add(StandardArgumentDefinitions.INPUT_LONG_NAME, inputBAM)
                .add(StandardArgumentDefinitions.OUTPUT_LONG_NAME, tempVcf.getPath())
                .add(StandardArgumentDefinitions.REFERENCE_LONG_NAME, b37_reference_20_21)
                .add("inference-batch-size", "2")
                .add("transfer-batch-size", "2")
                .add("tensor-type", tt.name())
                .add(StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false");

        runCommandLine(argsBuilder);
        assertInfoFieldsAreClose(tempVcf, expectedVcf, GATKVCFConstants.CNN_2D_KEY);
    }

    /**
     * Run the 2D Model on a small test VCF.
     */
    @Test(groups = {"python"})
    public void testInferenceArchitecture2d() {
        Utils.resetRandomGenerator();
        final boolean newExpectations = false;
        TensorType tt = TensorType.read_tensor;
        final String expectedVCFName = largeFileTestDir + "VQSR/expected/cnn_2d_chr20_subset_expected.vcf";
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        argsBuilder.add(StandardArgumentDefinitions.VARIANT_LONG_NAME, inputVCF)
                .add(StandardArgumentDefinitions.INPUT_LONG_NAME, inputBAM)
                .add(StandardArgumentDefinitions.REFERENCE_LONG_NAME, b37_reference_20_21)
                .add("architecture", architecture2D)
                .add("tensor-type", tt.name())
                .add("inference-batch-size", "8")
                .add("transfer-batch-size", "8")
                .add(StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false");

        if (newExpectations) {
            argsBuilder.add(StandardArgumentDefinitions.OUTPUT_LONG_NAME, expectedVCFName);
            runCommandLine(argsBuilder);
        } else {
            final File tempVcf = createTempFile("tester", ".vcf");
            final File expectedVcf = new File(expectedVCFName);
            argsBuilder.add(StandardArgumentDefinitions.OUTPUT_LONG_NAME, tempVcf.getPath());
            runCommandLine(argsBuilder);
            assertInfoFieldsAreClose(tempVcf, expectedVcf, GATKVCFConstants.CNN_2D_KEY);
        }
    }

    @Test(groups = {"python"})
    public void testInferenceWeights2d() {
        Utils.resetRandomGenerator();
        TensorType tt = TensorType.read_tensor;
        final File tempVcf = createTempFile("tester", ".vcf");
        final File expectedVcf = new File(largeFileTestDir + "VQSR/expected/cnn_2d_chr20_subset_expected.vcf");
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        argsBuilder.add(StandardArgumentDefinitions.VARIANT_LONG_NAME, inputVCF)
                .add(StandardArgumentDefinitions.INPUT_LONG_NAME, inputBAM)
                .add(StandardArgumentDefinitions.OUTPUT_LONG_NAME, tempVcf.getPath())
                .add(StandardArgumentDefinitions.REFERENCE_LONG_NAME, b37_reference_20_21)
                .add("weights", weights2D)
                .add("inference-batch-size", "4")
                .add("transfer-batch-size", "4")
                .add("tensor-type", tt.name())
                .add(StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false");

        runCommandLine(argsBuilder);
        assertInfoFieldsAreClose(tempVcf, expectedVcf, GATKVCFConstants.CNN_2D_KEY);
    }

    @Test(groups = {"python"})
    public void testInferenceArchitectureAndWeights2d() {
        Utils.resetRandomGenerator();
        TensorType tt = TensorType.read_tensor;
        final File tempVcf = createTempFile("tester", ".vcf");
        final File expectedVcf = new File(largeFileTestDir + "VQSR/expected/cnn_2d_chr20_subset_expected.vcf");
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        argsBuilder.add(StandardArgumentDefinitions.VARIANT_LONG_NAME, inputVCF)
                .add(StandardArgumentDefinitions.INPUT_LONG_NAME, inputBAM)
                .add(StandardArgumentDefinitions.OUTPUT_LONG_NAME, tempVcf.getPath())
                .add(StandardArgumentDefinitions.REFERENCE_LONG_NAME, b37_reference_20_21)
                .add("weights", weights2D)
                .add("architecture", architecture2D)
                .add("inference-batch-size", "4")
                .add("transfer-batch-size", "4")
                .add("tensor-type", tt.name())
                .add(StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false");

        runCommandLine(argsBuilder);
        assertInfoFieldsAreClose(tempVcf, expectedVcf, GATKVCFConstants.CNN_2D_KEY);
    }

    private void assertInfoFieldsAreClose(File actualVcf, File expectedVcf, String infoKey){
        Iterator<VariantContext> expectedVi = VariantContextTestUtils.streamVcf(expectedVcf).collect(Collectors.toList()).iterator();
        Iterator<VariantContext> actualVi = VariantContextTestUtils.streamVcf(actualVcf).collect(Collectors.toList()).iterator();
        while (expectedVi.hasNext() && actualVi.hasNext()) {
            VariantContext expectedVc = expectedVi.next();
            VariantContext actualVc = actualVi.next();
            double expectedScore = expectedVc.getAttributeAsDouble(infoKey, 0.0); // Different defaults trigger failures on missing scores
            double actualScore = actualVc.getAttributeAsDouble(infoKey, EPSILON+1.0);
            double diff = Math.abs(expectedScore-actualScore);
            Assert.assertTrue(diff < EPSILON);
            VariantContextTestUtils.assertVariantContextsAreEqual(actualVc, expectedVc, Collections.singletonList(infoKey), Collections.emptyList());
        }
        Assert.assertTrue(!expectedVi.hasNext() && !actualVi.hasNext());
    }



}
