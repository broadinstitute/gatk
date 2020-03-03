package org.broadinstitute.hellbender.tools.walkers.validation;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.engine.AbstractConcordanceWalker;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.nio.file.Path;

public class EvaluateInfoFieldConcordanceIntegrationTest extends CommandLineProgramTest {
    final double epsilon = 1e-3;

    @Test(dataProvider= "infoConcordanceDataProvider")
    public void testInfoConcordanceFromProvider(String inputVcf1, String inputVcf2, String evalKey, String truthKey,
                                                double snpMean, double snpSTD,
                                                double indelMean, double indelSTD) throws Exception {
        final Path summary = createTempPath("summary", ".txt");
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        argsBuilder.add(AbstractConcordanceWalker.EVAL_VARIANTS_SHORT_NAME, inputVcf1)
                .add(AbstractConcordanceWalker.TRUTH_VARIANTS_LONG_NAME, inputVcf2)
                .add("eval-info-key", evalKey)
                .add("truth-info-key", truthKey)
                .add(EvaluateInfoFieldConcordance.SUMMARY_LONG_NAME, summary.toString());
        runCommandLine(argsBuilder);

        try(InfoConcordanceRecord.InfoConcordanceReader
                    reader = new InfoConcordanceRecord.InfoConcordanceReader(summary)) {
            InfoConcordanceRecord snpRecord = reader.readRecord();
            InfoConcordanceRecord indelRecord = reader.readRecord();

            Assert.assertEquals(snpRecord.getVariantType(), VariantContext.Type.SNP);
            Assert.assertEquals(indelRecord.getVariantType(), VariantContext.Type.INDEL);

            Assert.assertEquals(snpRecord.getMean(), snpMean, epsilon);
            Assert.assertEquals(snpRecord.getStd(), snpSTD, epsilon);
            Assert.assertEquals(indelRecord.getMean(), indelMean, epsilon);
            Assert.assertEquals(indelRecord.getStd(), indelSTD, epsilon);
        }
    }

    @DataProvider
    public Object[][] infoConcordanceDataProvider() {
        return new Object [][]{
                new Object[]{
                        largeFileTestDir + "VQSR/expected/chr20_tiny_tf_python_gpu2.vcf",
                        largeFileTestDir + "VQSR/expected/chr20_tiny_tf_python_gpu2.vcf",
                        GATKVCFConstants.CNN_2D_KEY,
                        "NOVA_HISEQ_MIX_SMALL",
                        0.108878, 0.229415, 0.067024, 0.142705 // numbers verified by manual inspection

                },
                new Object[]{
                        largeFileTestDir + "VQSR/expected/chr20_tiny_tf_python_cpu.vcf",
                        largeFileTestDir + "VQSR/expected/chr20_tiny_th_python_gpu.vcf",
                        GATKVCFConstants.CNN_1D_KEY,
                        "NOVA_HISEQ_MIX_1D_RAB",
                        0.000256, 0.000136, 0.000240, 0.000153 // numbers verified by manual inspection
                }
        };
    }
}
