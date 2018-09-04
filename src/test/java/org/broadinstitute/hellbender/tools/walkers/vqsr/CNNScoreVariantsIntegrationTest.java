package org.broadinstitute.hellbender.tools.walkers.vqsr;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.annotations.Test;

import java.io.IOException;
import java.util.Arrays;

/**
 * Integration tests for {@link CNNScoreVariants}.
 * Created by sam on 1/8/18.
 */
public class CNNScoreVariantsIntegrationTest extends CommandLineProgramTest {
    private static final String architecture1D = largeFileTestDir + "VQSR/1d_cnn_mix_train_full_bn.json";
    private static final String architecture2D = largeFileTestDir + "VQSR/small_2d.json";
    private static final String inputVCF = largeFileTestDir + "VQSR/recalibrated_chr20_start.vcf";

    /**
     * Run the tool on a small test VCF.
     */
    @Test(groups = {"python"})
    public void testInference() throws IOException{
        final boolean newExpectations = false;
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        argsBuilder.addArgument(StandardArgumentDefinitions.VARIANT_LONG_NAME, inputVCF)
                .addArgument(StandardArgumentDefinitions.REFERENCE_LONG_NAME, b37_reference_20_21)
                .addArgument("architecture", architecture1D)
                .addArgument(StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false");

        if (newExpectations) {
            argsBuilder.addArgument(StandardArgumentDefinitions.OUTPUT_LONG_NAME, largeFileTestDir + "VQSR/expected/cnn_1d_chr20_subset_expected.vcf");
            runCommandLine(argsBuilder);
        } else {
            argsBuilder.addArgument(StandardArgumentDefinitions.OUTPUT_LONG_NAME, "%s");
            final IntegrationTestSpec spec = new IntegrationTestSpec(argsBuilder.toString(),
                    Arrays.asList(largeFileTestDir + "VQSR/expected/cnn_1d_chr20_subset_expected.vcf"));
            spec.executeTest("testInference", this);
        }
    }

    @Test(groups = {"python"})
    public void testInferenceWithWeightOverride() throws IOException{
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        argsBuilder.addArgument(StandardArgumentDefinitions.VARIANT_LONG_NAME, inputVCF)
                .addArgument(StandardArgumentDefinitions.OUTPUT_LONG_NAME, "%s")
                .addArgument(StandardArgumentDefinitions.REFERENCE_LONG_NAME, b37_reference_20_21)
                .addArgument("architecture", architecture1D)
                .addArgument("weights", architecture1D.replace(".json", ".hd5"))
                .addArgument(StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false");

        final IntegrationTestSpec spec = new IntegrationTestSpec(argsBuilder.toString(),
                Arrays.asList(largeFileTestDir + "VQSR/expected/cnn_1d_chr20_subset_expected.vcf"));
        spec.executeTest("testInference", this);
    }

    @Test(groups = {"python"})
    public void testInferenceWithWeightsOnly() throws IOException{
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        argsBuilder.addArgument(StandardArgumentDefinitions.VARIANT_LONG_NAME, inputVCF)
                .addArgument(StandardArgumentDefinitions.OUTPUT_LONG_NAME, "%s")
                .addArgument(StandardArgumentDefinitions.REFERENCE_LONG_NAME, b37_reference_20_21)
                .addArgument("weights", architecture1D.replace(".json", ".hd5"))
                .addArgument(StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false");

        final IntegrationTestSpec spec = new IntegrationTestSpec(argsBuilder.toString(),
                Arrays.asList(largeFileTestDir + "VQSR/expected/cnn_1d_chr20_subset_expected.vcf"));
        spec.executeTest("testInference", this);
    }

    @Test(groups = {"python"}, enabled = true)
    public void testInferenceResourceModel() throws IOException{
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        argsBuilder.addArgument(StandardArgumentDefinitions.VARIANT_LONG_NAME, inputVCF)
                .addArgument(StandardArgumentDefinitions.OUTPUT_LONG_NAME, "%s")
                .addArgument(StandardArgumentDefinitions.REFERENCE_LONG_NAME, b37_reference_20_21)
                .addArgument(StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false");

        final IntegrationTestSpec spec = new IntegrationTestSpec(argsBuilder.toString(),
                Arrays.asList(largeFileTestDir + "VQSR/expected/cnn_1d_chr20_subset_expected.vcf"));
        spec.executeTest("testInference", this);
    }

    @Test(groups = {"python"})
    public void testSmallBatchInference()throws IOException {
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        argsBuilder.addArgument(StandardArgumentDefinitions.VARIANT_LONG_NAME, inputVCF)
                .addArgument(StandardArgumentDefinitions.OUTPUT_LONG_NAME, "%s")
                .addArgument(StandardArgumentDefinitions.REFERENCE_LONG_NAME, b37_reference_20_21)
                .addArgument("architecture", architecture1D)
                .addArgument("inference-batch-size", "8")
                .addArgument("transfer-batch-size", "16")
                .addArgument(StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false");

        final IntegrationTestSpec spec = new IntegrationTestSpec(argsBuilder.toString(),
                Arrays.asList(largeFileTestDir + "VQSR/expected/cnn_1d_chr20_subset_expected.vcf"));
        spec.executeTest("testInference", this);
    }

    @Test(groups = {"python"})
    public void testOnContigEdge() throws IOException{
        final String edgeVcf = toolsTestDir + "walkers/VQSR/variantNearContigEdge.vcf";
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        argsBuilder.addArgument(StandardArgumentDefinitions.VARIANT_LONG_NAME, edgeVcf)
                .addArgument(StandardArgumentDefinitions.REFERENCE_LONG_NAME, hg19MiniReference)
                .addArgument("architecture", architecture1D)
                .addArgument(StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false");

        argsBuilder.addArgument(StandardArgumentDefinitions.OUTPUT_LONG_NAME, largeFileTestDir + "VQSR/expected/chrM.vcf");
        runCommandLine(argsBuilder);
    }

    /**
     * Run the 2D Model on a small test VCF.
     */
    @Test(groups = {"python"})
    public void testInference2d() throws IOException{
        final boolean newExpectations = false;
        TensorType tt = TensorType.read_tensor;
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        argsBuilder.addArgument(StandardArgumentDefinitions.VARIANT_LONG_NAME, inputVCF)
                .addArgument(StandardArgumentDefinitions.INPUT_LONG_NAME, largeFileTestDir + "VQSR/g94982_chr20_1m_10m_bamout.bam")
                .addArgument(StandardArgumentDefinitions.REFERENCE_LONG_NAME, b37_reference_20_21)
                .addArgument("architecture", architecture2D)
                .addArgument("inference-batch-size", "1")
                .addArgument("transfer-batch-size", "1")
                .addArgument("tensor-type", tt.name())
                .addArgument(StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false");

        if (newExpectations) {
            argsBuilder.addArgument(StandardArgumentDefinitions.OUTPUT_LONG_NAME, largeFileTestDir + "VQSR/expected/cnn_2d_chr20_subset_expected.vcf");
            runCommandLine(argsBuilder);
        } else {
            argsBuilder.addArgument(StandardArgumentDefinitions.OUTPUT_LONG_NAME, "%s");
            final IntegrationTestSpec spec = new IntegrationTestSpec(argsBuilder.toString(),
                    Arrays.asList(largeFileTestDir + "VQSR/expected/cnn_2d_chr20_subset_expected.vcf"));
            spec.executeTest("testInference2d", this);
        }

    }

    /**
     * Run the 2D Model on a small test VCF with the resource loaded architecture.
     */
    @Test(groups = {"python"}, enabled = true)
    public void testInference2dResourceModel() throws IOException{
        TensorType tt = TensorType.read_tensor;
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        argsBuilder.addArgument(StandardArgumentDefinitions.VARIANT_LONG_NAME, inputVCF)
                .addArgument(StandardArgumentDefinitions.INPUT_LONG_NAME, largeFileTestDir + "VQSR/g94982_chr20_1m_10m_bamout.bam")
                .addArgument(StandardArgumentDefinitions.OUTPUT_LONG_NAME, "%s")
                .addArgument(StandardArgumentDefinitions.REFERENCE_LONG_NAME, b37_reference_20_21)
                .addArgument("inference-batch-size", "1")
                .addArgument("transfer-batch-size", "1")
                .addArgument("tensor-type", tt.name())
                .addArgument(StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false");

        final IntegrationTestSpec spec = new IntegrationTestSpec(argsBuilder.toString(),
                Arrays.asList(largeFileTestDir + "VQSR/expected/cnn_2d_chr20_subset_expected.vcf"));
        spec.executeTest("testInference2d", this);

    }
}
