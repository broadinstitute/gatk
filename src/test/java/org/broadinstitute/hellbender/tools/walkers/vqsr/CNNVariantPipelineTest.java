package org.broadinstitute.hellbender.tools.walkers.vqsr;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.Main;
import org.testng.annotations.Test;

import java.io.File;


public class CNNVariantPipelineTest extends GATKBaseTest {
    final private static String inputVCF = largeFileTestDir + "VQSR/g94982_b37_chr20_1m_10m.vcf";
    final private static String truthVCF = largeFileTestDir + "VQSR/giab_chr20_1m_10m.vcf.gz";
    final private static String truthBED = largeFileTestDir + "VQSR/giab_na12878_confident_chr20_1m_10m.bed";

    private static File readTensorDir;
    private static File referenceTensorDir;

    @Test(groups = {"python"})
    public static void makeTempDirectories() {
        readTensorDir = createTempDir("readTensorDir");
        referenceTensorDir = createTempDir("referenceTensorDir");
    }

    @Test(groups = {"python"}, dependsOnMethods = {"makeTempDirectories"})
    public void testGenerateReferenceTensors() {
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

    @Test(groups = {"python"}, dependsOnMethods = {"makeTempDirectories"})
    public void testGenerateReadTensors() {
        final String bamFile = largeFileTestDir + "VQSR/g94982_chr20_1m_10m_bamout.bam";
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

    @Test(groups = {"python"}, dependsOnMethods = {"testGenerateReferenceTensors"})
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

    @Test(groups = {"python"}, dependsOnMethods = {"testGenerateReadTensors"})
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
        final String trancheVCF = largeFileTestDir + "VQSR/g94982_b37_chr20_1m_10m.vcf.gz";
        final String snpTruthVCF = largeFileTestDir + "VQSR/giab_chr20_1m_10m.vcf.gz";
        final String indelTruthVCF = largeFileTestDir + "VQSR/giab_chr20_1m_10m.vcf.gz";
        final File outputVCF = createTempFile("variant_tranches_output", "vcf");
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("FilterVariantTranches")
                .addArgument(StandardArgumentDefinitions.VARIANT_LONG_NAME, trancheVCF)
                .addArgument(StandardArgumentDefinitions.OUTPUT_LONG_NAME, outputVCF.getAbsolutePath())
                .addArgument("resource", snpTruthVCF)
                .addArgument("resource", indelTruthVCF)
                .addArgument("snp-tranche", "99.9")
                .addArgument("indel-tranche", "99.0")
                .addArgument("info-key", "VQSLOD");

        new Main().instanceMain(args.getArgsArray());
    }
}
