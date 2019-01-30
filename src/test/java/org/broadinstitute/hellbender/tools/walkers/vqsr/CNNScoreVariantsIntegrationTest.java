package org.broadinstitute.hellbender.tools.walkers.vqsr;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.broadinstitute.hellbender.utils.Utils;
import org.testng.annotations.Test;

import java.io.IOException;
import java.util.Arrays;

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

    /**
     * Run the tool on a small test VCF.
     */
    @Test(groups = {"python"})
    public void testAllDefaultArgs() throws IOException {
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
    public void testInferenceArchitecture() throws IOException {
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
    public void testInferenceWeights() throws IOException {
        final boolean newExpectations = false;
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        argsBuilder.addArgument(StandardArgumentDefinitions.VARIANT_LONG_NAME, inputVCF)
                .addArgument(StandardArgumentDefinitions.REFERENCE_LONG_NAME, b37_reference_20_21)
                .addArgument("weights", weights1D)
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
    public void testInferenceArchitectureAndWeights() throws IOException {
        final boolean newExpectations = false;
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        argsBuilder.addArgument(StandardArgumentDefinitions.VARIANT_LONG_NAME, inputVCF)
                .addArgument(StandardArgumentDefinitions.REFERENCE_LONG_NAME, b37_reference_20_21)
                .addArgument("weights", weights1D)
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
    public void testInferenceWithIntervals() throws IOException {
        final boolean newExpectations = false;
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        argsBuilder.addArgument(StandardArgumentDefinitions.VARIANT_LONG_NAME, bigInputVCF)
                .addArgument(StandardArgumentDefinitions.REFERENCE_LONG_NAME, b37_reference_20_21)
                .addArgument(StandardArgumentDefinitions.INTERVALS_LONG_NAME, inputIntervals)
                .addArgument(StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false");

        if (newExpectations) {
            argsBuilder.addArgument(StandardArgumentDefinitions.OUTPUT_LONG_NAME, largeFileTestDir + "VQSR/expected/cnn_1d_contig20_1m_10m_expected.vcf");
            runCommandLine(argsBuilder);
        } else {
            argsBuilder.addArgument(StandardArgumentDefinitions.OUTPUT_LONG_NAME, "%s");
            final IntegrationTestSpec spec = new IntegrationTestSpec(argsBuilder.toString(),
                    Arrays.asList(largeFileTestDir + "VQSR/expected/cnn_1d_contig20_1m_10m_expected.vcf"));
            spec.executeTest("testInference", this);
        }
    }

    @Test(groups = {"python"})
    public void testSmallBatchInference() throws IOException {
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        argsBuilder.addArgument(StandardArgumentDefinitions.VARIANT_LONG_NAME, inputVCF)
                .addArgument(StandardArgumentDefinitions.OUTPUT_LONG_NAME, "%s")
                .addArgument(StandardArgumentDefinitions.REFERENCE_LONG_NAME, b37_reference_20_21)
                .addArgument("inference-batch-size", "8")
                .addArgument("transfer-batch-size", "16")
                .addArgument(StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false");

        final IntegrationTestSpec spec = new IntegrationTestSpec(argsBuilder.toString(),
                Arrays.asList(largeFileTestDir + "VQSR/expected/cnn_1d_chr20_subset_expected.vcf"));
        spec.executeTest("testInference", this);
    }

    @Test(groups = {"python"})
    public void testOnContigEdge() throws IOException {
        final String edgeVcf = toolsTestDir + "walkers/VQSR/variantNearContigEdge.vcf";
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        argsBuilder.addArgument(StandardArgumentDefinitions.VARIANT_LONG_NAME, edgeVcf)
                .addArgument(StandardArgumentDefinitions.OUTPUT_LONG_NAME, "%s")
                .addArgument(StandardArgumentDefinitions.REFERENCE_LONG_NAME, hg19MiniReference)
                .addArgument(StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false");

        final IntegrationTestSpec spec = new IntegrationTestSpec(argsBuilder.toString(),
                Arrays.asList(largeFileTestDir + "VQSR/expected/chrM.vcf"));
        spec.executeTest("testContigOnEdge", this);
    }

    /**
     * Run the 2D Model on a small test VCF with the resource loaded weights and architecture.
     */
    @Test(groups = {"python"})
    public void testInference2dResourceModel() throws IOException {
        // We reset the random number generator at the beginning of each test so that the random down-sampling of reads
        // by the reservoir down-sampler does not cause slightly different scores.
        Utils.getRandomGenerator().setSeed(1234);
        TensorType tt = TensorType.read_tensor;
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        argsBuilder.addArgument(StandardArgumentDefinitions.VARIANT_LONG_NAME, inputVCF)
                .addArgument(StandardArgumentDefinitions.INPUT_LONG_NAME, inputBAM)
                .addArgument(StandardArgumentDefinitions.OUTPUT_LONG_NAME, "%s")
                .addArgument(StandardArgumentDefinitions.REFERENCE_LONG_NAME, b37_reference_20_21)
                .addArgument("inference-batch-size", "2")
                .addArgument("transfer-batch-size", "2")
                .addArgument("tensor-type", tt.name())
                .addArgument(StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false");

        final IntegrationTestSpec spec = new IntegrationTestSpec(argsBuilder.toString(),
                Arrays.asList(largeFileTestDir + "VQSR/expected/cnn_2d_chr20_subset_expected.vcf"));
        spec.executeTest("testInference2d", this);

    }

    /**
     * Run the 2D Model on a small test VCF.
     */
    @Test(groups = {"python"})
    public void testInferenceArchitecture2d() throws IOException {
        Utils.getRandomGenerator().setSeed(1234);
        final boolean newExpectations = false;
        TensorType tt = TensorType.read_tensor;
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        argsBuilder.addArgument(StandardArgumentDefinitions.VARIANT_LONG_NAME, inputVCF)
                .addArgument(StandardArgumentDefinitions.INPUT_LONG_NAME, inputBAM)
                .addArgument(StandardArgumentDefinitions.REFERENCE_LONG_NAME, b37_reference_20_21)
                .addArgument("architecture", architecture2D)
                .addArgument("inference-batch-size", "4")
                .addArgument("transfer-batch-size", "4")
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

    @Test(groups = {"python"})
    public void testInferenceWeights2d() throws IOException {
        Utils.getRandomGenerator().setSeed(1234);
        TensorType tt = TensorType.read_tensor;
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        argsBuilder.addArgument(StandardArgumentDefinitions.VARIANT_LONG_NAME, inputVCF)
                .addArgument(StandardArgumentDefinitions.INPUT_LONG_NAME, inputBAM)
                .addArgument(StandardArgumentDefinitions.OUTPUT_LONG_NAME, "%s")
                .addArgument(StandardArgumentDefinitions.REFERENCE_LONG_NAME, b37_reference_20_21)
                .addArgument("weights", weights2D)
                .addArgument("inference-batch-size", "4")
                .addArgument("transfer-batch-size", "4")
                .addArgument("tensor-type", tt.name())
                .addArgument(StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false");


        final IntegrationTestSpec spec = new IntegrationTestSpec(argsBuilder.toString(),
                    Arrays.asList(largeFileTestDir + "VQSR/expected/cnn_2d_chr20_subset_expected.vcf"));
        spec.executeTest("testInference2d", this);
    }

    @Test(groups = {"python"})
    public void testInferenceArchitectureAndWeights2d() throws IOException {
        Utils.getRandomGenerator().setSeed(1234);
        TensorType tt = TensorType.read_tensor;
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        argsBuilder.addArgument(StandardArgumentDefinitions.VARIANT_LONG_NAME, inputVCF)
                .addArgument(StandardArgumentDefinitions.INPUT_LONG_NAME, inputBAM)
                .addArgument(StandardArgumentDefinitions.OUTPUT_LONG_NAME, "%s")
                .addArgument(StandardArgumentDefinitions.REFERENCE_LONG_NAME, b37_reference_20_21)
                .addArgument("weights", weights2D)
                .addArgument("architecture", architecture2D)
                .addArgument("inference-batch-size", "4")
                .addArgument("transfer-batch-size", "4")
                .addArgument("tensor-type", tt.name())
                .addArgument(StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false");


        final IntegrationTestSpec spec = new IntegrationTestSpec(argsBuilder.toString(),
                Arrays.asList(largeFileTestDir + "VQSR/expected/cnn_2d_chr20_subset_expected.vcf"));
        spec.executeTest("testInference2d", this);
    }

}
