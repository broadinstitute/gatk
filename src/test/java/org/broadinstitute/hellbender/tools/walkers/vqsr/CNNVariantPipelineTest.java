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
        args.addRaw("CNNVariantWriteTensors")
                .add(StandardArgumentDefinitions.VARIANT_LONG_NAME, inputVCF)
                .add(StandardArgumentDefinitions.REFERENCE_LONG_NAME, b37_reference_20_21)
                .add("truth-vcf", truthVCF)
                .add("truth-bed", truthBED)
                .add("tensor-type", TensorType.reference.name())
                .add("output-tensor-dir", referenceTensorDir.toString());

        new Main().instanceMain(args.getArgsArray());
    }

    @Test(groups = {"python"}, dependsOnMethods = {"makeTempDirectories"})
    public void testGenerateReadTensors() {
        final String bamFile = largeFileTestDir + "VQSR/g94982_b37_chr20_1m_8m_bamout.bam";
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addRaw("CNNVariantWriteTensors")
                .add(StandardArgumentDefinitions.VARIANT_LONG_NAME, inputVCF)
                .add(StandardArgumentDefinitions.REFERENCE_LONG_NAME, b37_reference_20_21)
                .add("truth-vcf", truthVCF)
                .add("truth-bed", truthBED)
                .add("bam-file", bamFile)
                .add("max-tensors", "4000")
                .add("tensor-type", TensorType.read_tensor.name())
                .add("output-tensor-dir", readTensorDir.toString())
                .add("channels-last", "true");

        new Main().instanceMain(args.getArgsArray());
    }

    @Test(groups = {"python"}, dependsOnMethods = {"testGenerateReferenceTensors"})
    public void testTrainingReferenceModel() {
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addRaw("CNNVariantTrain")
                .add("input-tensor-dir", referenceTensorDir.toString()+"/")
                .add("tensor-type", TensorType.reference.name())
                .add("epochs", "1")
                .add("training-steps", "30")
                .add("model-name", "test_reference_model")
                .add("output-dir", referenceTensorDir.toString()+"/");

        new Main().instanceMain(args.getArgsArray());
    }

    @Test(groups = {"python"}, dependsOnMethods = {"testGenerateReadTensors"})
    public void testTrainingReadModel() {
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addRaw("CNNVariantTrain")
                .add("input-tensor-dir", readTensorDir.toString()+"/")
                .add("tensor-type", TensorType.read_tensor.name())
                .add("epochs", "1")
                .add("training-steps", "5")
                .add("validation-steps", "2")
                .add("model-name", "test_read_tensor_model")
                .add("output-dir", readTensorDir.toString()+"/")
                .add("channels-last", "true");

        new Main().instanceMain(args.getArgsArray());
    }

    @Test(groups = {"python"})
    public void testTranches() {
        final String trancheVCF = largeFileTestDir + "VQSR/g94982_b37_chr20_1m_10m.vcf.gz";
        final String snpTruthVCF = largeFileTestDir + "VQSR/giab_chr20_1m_10m.vcf.gz";
        final String indelTruthVCF = largeFileTestDir + "VQSR/giab_chr20_1m_10m.vcf.gz";
        final File outputVCF = createTempFile("variant_tranches_output", "vcf");
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addRaw("FilterVariantTranches")
                .add(StandardArgumentDefinitions.VARIANT_LONG_NAME, trancheVCF)
                .add(StandardArgumentDefinitions.OUTPUT_LONG_NAME, outputVCF.getAbsolutePath())
                .add(StandardArgumentDefinitions.RESOURCE_LONG_NAME, snpTruthVCF)
                .add(StandardArgumentDefinitions.RESOURCE_LONG_NAME, indelTruthVCF)
                .add("snp-tranche", "99.9")
                .add("indel-tranche", "99.0")
                .add("info-key", "VQSLOD");

        new Main().instanceMain(args.getArgsArray());
    }
}
