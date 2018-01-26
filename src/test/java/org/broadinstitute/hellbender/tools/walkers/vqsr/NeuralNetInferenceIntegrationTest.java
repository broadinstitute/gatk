package org.broadinstitute.hellbender.tools.walkers.vqsr;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.test.IntegrationTestSpec;
import org.broadinstitute.hellbender.utils.variant.VcfUtils;
import org.testng.annotations.Test;

import java.io.IOException;
import java.util.Arrays;

/**
 * Integration tests for {@link NeuralNetInference}.
 * Created by sam on 1/8/18.
 */
public class NeuralNetInferenceIntegrationTest extends CommandLineProgramTest {
    private static String architectureHD5 = packageMainResourcesDir + "tools/walkers/vqsr/cnn_1d_annotations.hd5";

    /**
     * Run the tool on a small test VCF.
     */
    @Test(groups = {"python"})
    public void testInference() throws IOException{
        final String inputVCF = largeFileTestDir + "VQSR/recalibrated_chr20_start.vcf";
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        argsBuilder.addArgument(StandardArgumentDefinitions.VARIANT_LONG_NAME, inputVCF)
                .addArgument(StandardArgumentDefinitions.OUTPUT_LONG_NAME, "%s")
                .addArgument(StandardArgumentDefinitions.REFERENCE_LONG_NAME, b37_reference_20_21)
                .addArgument("architecture", architectureHD5)
                .addArgument(StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false");

        final IntegrationTestSpec spec = new IntegrationTestSpec(argsBuilder.toString(),
                Arrays.asList(largeFileTestDir + "VQSR/expected/cnn_1d_chr20_subset_expected.vcf"));
        spec.executeTest("testInference", this);

    }

//    @Test(groups = {"python"})
//    public void generateTestExample() throws IOException{
//        final String inputVCF = largeFileTestDir + "VQSR/recalibrated_chr20_start.vcf";
//        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
//        argsBuilder.addArgument(StandardArgumentDefinitions.VARIANT_LONG_NAME, inputVCF)
//                .addArgument(StandardArgumentDefinitions.OUTPUT_LONG_NAME, "cnn_1d_chr20_subset_expected.vcf")
//                .addArgument(StandardArgumentDefinitions.REFERENCE_LONG_NAME, b37_reference_20_21)
//                .addArgument("architecture", architectureHD5)
//                .addArgument(StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false");
//
//        runCommandLine(argsBuilder);
//
//    }


    @Test(groups = {"python"})
    public void testSmallBatchInference()throws IOException {
        final String inputVCF = largeFileTestDir + "VQSR/recalibrated_chr20_start.vcf";
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        argsBuilder.addArgument(StandardArgumentDefinitions.VARIANT_LONG_NAME, inputVCF)
                .addArgument(StandardArgumentDefinitions.OUTPUT_LONG_NAME, "%s")
                .addArgument(StandardArgumentDefinitions.REFERENCE_LONG_NAME, b37_reference_20_21)
                .addArgument("architecture", architectureHD5)
                .addArgument("inference-batch-size", "8")
                .addArgument("transfer-batch-size", "16")
                .addArgument(StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false");

        final IntegrationTestSpec spec = new IntegrationTestSpec(argsBuilder.toString(),
                Arrays.asList(largeFileTestDir + "VQSR/expected/cnn_1d_chr20_subset_expected.vcf"));
        spec.executeTest("testInference", this);
    }
}
