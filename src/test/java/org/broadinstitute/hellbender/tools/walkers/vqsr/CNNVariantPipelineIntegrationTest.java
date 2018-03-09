package org.broadinstitute.hellbender.tools.walkers.vqsr;

import org.broadinstitute.hellbender.Main;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;

import static org.broadinstitute.hellbender.GATKBaseTest.b37_reference_20_21;
import static org.broadinstitute.hellbender.GATKBaseTest.largeFileTestDir;

public class CNNVariantPipelineIntegrationTest {
    final private static String inputVCF = largeFileTestDir + "VQSR/g94982_b37_chr20_1m_10m.vcf";
    final private static String truthVCF = largeFileTestDir + "VQSR/giab_chr20_1m_10m.vcf.gz";
    final private static String truthBED = largeFileTestDir + "VQSR/giab_na12878_confident_chr20_1m_10m.bed";
    final private static String bamFile = largeFileTestDir + "VQSR/g94982_chr20_1m_10m_bamout.bam";

    final private static String trancheVCF = largeFileTestDir + "VQSR/g94982_b37_chr20_1m_10m.vcf.gz";
    final private static String snpTruthVCF = largeFileTestDir + "VQSR/giab_chr20_1m_10m.vcf.gz";
    final private static String indelTruthVCF = largeFileTestDir + "VQSR/giab_chr20_1m_10m.vcf.gz";
    final private static String outputVCF = largeFileTestDir + "VQSR/expected/variant_tranches_python_expected.vcf";

    private static Path readTensorDir;
    private static Path referenceTensorDir;


    @BeforeClass(alwaysRun = true)
    public static void makeTempDirectories() throws IOException {
        readTensorDir = Files.createTempDirectory("readTensorDir");
        referenceTensorDir = Files.createTempDirectory("referenceTensorDir");
        readTensorDir.toFile().deleteOnExit();
        referenceTensorDir.toFile().deleteOnExit();
    }

    @Test(groups = {"python"})
    public void generateReferenceTensors() {
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("CNNVariantWriteTensors")
                .addArgument(StandardArgumentDefinitions.VARIANT_LONG_NAME, inputVCF)
                .addArgument(StandardArgumentDefinitions.REFERENCE_LONG_NAME, b37_reference_20_21)
                .addArgument("truth-vcf", truthVCF)
                .addArgument("truth-bed", truthBED)
                .addArgument("tensor-type", TensorType.reference.name())
                .addArgument("output-tensor-dir", referenceTensorDir.toString());

        new Main().instanceMain(args.getArgsArray());
    }

    @Test(groups = {"python"})
    public void generateReadTensors() {
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("CNNVariantWriteTensors")
                .addArgument(StandardArgumentDefinitions.VARIANT_LONG_NAME, inputVCF)
                .addArgument(StandardArgumentDefinitions.REFERENCE_LONG_NAME, b37_reference_20_21)
                .addArgument("truth-vcf", truthVCF)
                .addArgument("truth-bed", truthBED)
                .addArgument("bam-file", bamFile)
                .addArgument("max-tensors", "4000")
                .addArgument("tensor-type", TensorType.read_tensor.name())
                .addArgument("output-tensor-dir", readTensorDir.toString())
                .addArgument("channels-last", "true");

        new Main().instanceMain(args.getArgsArray());
    }

    @Test(groups = {"python"}, dependsOnMethods = {"generateReferenceTensors"})
    public void testTrainingReferenceModel() {
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("CNNVariantTrain")
                .addArgument("input-tensor-dir", referenceTensorDir.toString()+"/")
                .addArgument("tensor-type", TensorType.reference.name())
                .addArgument("epochs", "1")
                .addArgument("training-steps", "30")
                .addArgument("model-name", "test_reference_model")
                .addArgument("output-dir", referenceTensorDir.toString()+"/");

        new Main().instanceMain(args.getArgsArray());
    }

    @Test(groups = {"python"}, dependsOnMethods = {"generateReadTensors"})
    public void testTrainingReadModel() {
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("CNNVariantTrain")
                .addArgument("input-tensor-dir", readTensorDir.toString()+"/")
                .addArgument("tensor-type", TensorType.read_tensor.name())
                .addArgument("epochs", "1")
                .addArgument("training-steps", "5")
                .addArgument("validation-steps", "2")
                .addArgument("model-name", "test_read_tensor_model")
                .addArgument("output-dir", readTensorDir.toString()+"/")
                .addArgument("channels-last", "true");

        new Main().instanceMain(args.getArgsArray());
    }

    @Test(groups = {"python"})
    public void testTranches() {
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("VariantTranchesFromInfoKey")
                .addArgument(StandardArgumentDefinitions.VARIANT_LONG_NAME, trancheVCF)
                .addArgument(StandardArgumentDefinitions.OUTPUT_LONG_NAME, outputVCF)
                .addArgument("snp-truth-vcf", snpTruthVCF)
                .addArgument("indel-truth-vcf", indelTruthVCF)
                .addArgument("tranche", "99.0")
                .addArgument("tranche", "95.0")
                .addArgument("max-sites", "2000")
                .addArgument("info-key", "VQSLOD");

        new Main().instanceMain(args.getArgsArray());
    }
}
