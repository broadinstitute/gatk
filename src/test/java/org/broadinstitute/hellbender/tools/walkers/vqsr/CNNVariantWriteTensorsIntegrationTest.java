package org.broadinstitute.hellbender.tools.walkers.vqsr;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.test.IntegrationTestSpec;
import org.testng.annotations.Test;

import java.io.IOException;
import java.util.Arrays;

/**
 * Created by sam on 2/7/18.
 */
public class CNNVariantWriteTensorsIntegrationTest extends CommandLineProgramTest {
    final private static String inputVCF = largeFileTestDir + "VQSR/g94982_b37_chr20_1m_10m.vcf";
    final private static String truthVCF = largeFileTestDir + "VQSR/giab_chr20_1m_10m.vcf.gz";
    final private static String truthBED = largeFileTestDir + "VQSR/giab_na12878_confident_chr20_1m_10m.bed";
    final private static String bamFile = largeFileTestDir + "CEUTrio.HiSeq.WGS.b37.NA12878.20.21.bam";


    @Test(groups = {"python"})
    public void generateReferenceTensors() throws IOException{
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        argsBuilder.addArgument(StandardArgumentDefinitions.VARIANT_LONG_NAME, inputVCF)
                .addArgument(StandardArgumentDefinitions.REFERENCE_LONG_NAME, b37_reference_20_21)
                .addArgument("truth-vcf", truthVCF)
                .addArgument("truth-bed", truthBED)
                .addArgument("tensor-name", TensorMapEnum.reference.name())
                .addArgument("data-dir", largeFileTestDir + "VQSR/reference_tensors/");

        runCommandLine(argsBuilder);
    }

    @Test(groups = {"python"})
    public void generateReadTensors() throws IOException{
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        argsBuilder.addArgument(StandardArgumentDefinitions.VARIANT_LONG_NAME, inputVCF)
                .addArgument(StandardArgumentDefinitions.REFERENCE_LONG_NAME, b37_reference_20_21)
                .addArgument("truth-vcf", truthVCF)
                .addArgument("truth-bed", truthBED)
                .addArgument("bam-file", bamFile)
                .addArgument("max-tensors ", "1000")
                .addArgument("tensor-name", TensorMapEnum.reference.name())
                .addArgument("data-dir", largeFileTestDir + "VQSR/read_tensors/");

        runCommandLine(argsBuilder);
    }


}

