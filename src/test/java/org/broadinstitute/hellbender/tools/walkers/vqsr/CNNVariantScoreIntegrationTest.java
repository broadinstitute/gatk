package org.broadinstitute.hellbender.tools.walkers.vqsr;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.test.IntegrationTestSpec;
import org.testng.annotations.Test;

import java.io.IOException;
import java.util.Arrays;

/**
 * Integration tests for {@link CNNVariantScore}.
 * Created by sam on 1/8/18.
 */
public class CNNVariantScoreIntegrationTest extends CommandLineProgramTest {
    private static String architecture1D = packageMainResourcesDir + "tools/walkers/vqsr/1d_cnn_mix_train_full_bn.json";
    private static String architecture2D = packageMainResourcesDir + "tools/walkers/vqsr/2d_cnn_mix_train.json";

    private static final String inputVCF = largeFileTestDir + "VQSR/recalibrated_chr20_start.vcf";

    /**
     * Run the tool on a small test VCF.
     */
    @Test(groups = {"python"})
    public void testInference() throws IOException{
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        argsBuilder.addArgument(StandardArgumentDefinitions.VARIANT_LONG_NAME, inputVCF)
                .addArgument(StandardArgumentDefinitions.OUTPUT_LONG_NAME, "%s")
                .addArgument(StandardArgumentDefinitions.REFERENCE_LONG_NAME, b37_reference_20_21)
                .addArgument("architecture", architecture1D)
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
//                .addArgument(StandardArgumentDefinitions.OUTPUT_LONG_NAME, largeFileTestDir + "VQSR/expected/cnn_1d_chr20_subset_expected.vcf")
//                .addArgument(StandardArgumentDefinitions.REFERENCE_LONG_NAME, b37_reference_20_21)
//                .addArgument("architecture", architecture1D)
//                .addArgument(StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false");
//
//        runCommandLine(argsBuilder);
//    }


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


    /**
     * Run the 2D Model on a small test VCF.
     */
    @Test(groups = {"python"})
    public void testInference2d() throws IOException{
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        argsBuilder.addArgument(StandardArgumentDefinitions.VARIANT_LONG_NAME, inputVCF)
                .addArgument(StandardArgumentDefinitions.INPUT_LONG_NAME, largeFileTestDir + "CEUTrio.HiSeq.WGS.b37.NA12878.20.21.bam")
                .addArgument(StandardArgumentDefinitions.OUTPUT_LONG_NAME, "%s")
                .addArgument(StandardArgumentDefinitions.REFERENCE_LONG_NAME, b37_reference_20_21)
                .addArgument("architecture", architecture2D)
                .addArgument("inference-batch-size", "1")
                .addArgument("transfer-batch-size", "1")
                .addArgument("use-reads", "true")
                .addArgument(StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false");

        //runCommandLine(argsBuilder);
        final IntegrationTestSpec spec = new IntegrationTestSpec(argsBuilder.toString(),
                Arrays.asList(largeFileTestDir + "VQSR/expected/cnn_2d_chr20_subset_expected.vcf"));
        spec.executeTest("testInference2d", this);

    }
}
