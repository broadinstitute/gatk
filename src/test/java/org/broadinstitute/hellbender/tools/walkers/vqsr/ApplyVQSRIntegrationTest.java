package org.broadinstitute.hellbender.tools.walkers.vqsr;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.IntegrationTestSpec;
import org.testng.annotations.Test;

import java.io.IOException;
import java.util.Arrays;

public class ApplyVQSRIntegrationTest extends CommandLineProgramTest{

    @Test
    public void testMissingTranches() throws IOException {
        IntegrationTestSpec spec = new IntegrationTestSpec(
                " -mode BOTH" +
                        " --excludeFiltered -ts_filter_level 90.0" +
                        " -variant " + getToolTestDataDir() + "VQSR.mixedTest.input" +
                        " -vcfOut %s" +
                        " -tranchesFile " + getToolTestDataDir() + "VQSR.mixedTest.NO.tranches" + //this file exists
                        " -recalFile " + getToolTestDataDir() + "VQSR.mixedTest.recal",
                1,
                UserException.BadInput.class
        );
        spec.executeTest("testMissingTranches", this);
    }

    @Test
    public void testMissingTrancheFile() throws IOException {
        IntegrationTestSpec spec = new IntegrationTestSpec(
                " -mode BOTH" +
                        " --excludeFiltered -ts_filter_level 90.0" +
                        " -variant " + getToolTestDataDir() + "VQSR.mixedTest.input" +
                        " -vcfOut %s" +
                        " -tranchesFile " + getToolTestDataDir() + "VQSR.mixedTest.MISSING.tranches" +
                        " -recalFile " + getToolTestDataDir() + "VQSR.mixedTest.recal",
                1,
                UserException.CouldNotReadInputFile.class
        );
        spec.executeTest("testMissingTrancheFile", this);
    }

    @Test
    public void testBadOptionCombination() throws IOException {
        IntegrationTestSpec spec = new IntegrationTestSpec(
                " -mode BOTH" +
                        " --lodCutoff 0.0 -ts_filter_level 90.0" +
                        " -variant " + getToolTestDataDir() + "VQSR.mixedTest.input" +
                        " -vcfOut %s" +
                        " -tranchesFile " + getToolTestDataDir() + "VQSR.mixedTest.tranches" +
                        " -recalFile " + getToolTestDataDir() + "VQSR.mixedTest.recal",
                1,
                UserException.CommandLineException.class
        );
        spec.executeTest("testBadOptionCombination", this);
    }

    @Test
    public void testApplyRecalibrationSnpAndIndelTogetherNoIntervalsExcludeFiltered() throws IOException {
        IntegrationTestSpec spec = new IntegrationTestSpec(
                " -mode BOTH" +
                        " --excludeFiltered -ts_filter_level 90.0" +
                        " -variant " + getToolTestDataDir() + "VQSR.mixedTest.input" +
                        " -vcfOut %s" +
                        " -tranchesFile " + getToolTestDataDir() + "VQSR.mixedTest.tranches" +
                        " -recalFile " + getToolTestDataDir() + "VQSR.mixedTest.recal",
                Arrays.asList(getToolTestDataDir() + "expected.VQSR.mixedTest.noRanges.excludeFiltered_ts_filter_level90.vcf")
        );
        spec.executeTest("testApplyRecalibrationSnpAndIndelTogetherNoIntervalsExcludeFiltered", this);
    }

    @Test
    public void testApplyRecalibrationVariantNotInRecalFile() throws IOException {
        IntegrationTestSpec spec = new IntegrationTestSpec(
                " -mode INDEL" +
                        " -variant " + getToolTestDataDir() + "VQSR.mixedTest.unknownVariant.input" + //contains a variant not in the recal file
                        " -vcfOut %s" +
                        " -tranchesFile " + getToolTestDataDir() + "VQSR.mixedTest.tranches" +
                        " -recalFile " + getToolTestDataDir() + "VQSR.mixedTest.recal",
                1,
                UserException.BadInput.class
        );
        spec.executeTest("testApplyRecalibrationVariantNotInRecalFile", this);
    }

    @Test
    public void testApplyRecalibrationNoVQSLOD() throws IOException {
        IntegrationTestSpec spec = new IntegrationTestSpec(
                " -mode BOTH" +
                        " -variant " + getToolTestDataDir() + "VQSR.mixedTest.input" +
                        " -vcfOut %s" +
                        " -tranchesFile " + getToolTestDataDir() + "VQSR.mixedTest.tranches" +
                        " -recalFile " + getToolTestDataDir() + "VQSR.mixedTest.NO_VQSLOD.recal",
                1,
                UserException.BadInput.class
        );
        spec.executeTest("testApplyRecalibrationNoVQSLOD", this);
    }

    @Test
    public void testApplyRecalibrationBadVQSLOD() throws IOException {
        IntegrationTestSpec spec = new IntegrationTestSpec(
                " -mode BOTH" +
                        " -variant " + getToolTestDataDir() + "VQSR.mixedTest.input" +
                        " -vcfOut %s" +
                        " -tranchesFile " + getToolTestDataDir() + "VQSR.mixedTest.tranches" +
                        " -recalFile " + getToolTestDataDir() + "VQSR.mixedTest.BAD_VQSLOD.recal",
                1,
                UserException.BadInput.class
        );
        spec.executeTest("testApplyRecalibrationBadVQSLOD", this);
    }

    @Test
    public void testApplyRecalibrationSnpNoIntervals() throws IOException {
        IntegrationTestSpec spec = new IntegrationTestSpec(
                " -mode SNP" +
                        " -variant " + getToolTestDataDir() + "VQSR.mixedTest.input" +
                        " -vcfOut %s" +
                        " -tranchesFile " + getToolTestDataDir() + "VQSR.mixedTest.tranches" +
                        " -recalFile " + getToolTestDataDir() + "VQSR.mixedTest.recal",
                Arrays.asList(getToolTestDataDir() + "expected.VQSR.mixedTest.noRanges.modeSNP.vcf")
        );
        spec.executeTest("testApplyRecalibrationSnpNoIntervals", this);
    }

    @Test
    public void testApplyRecalibrationIndelNoIntervals() throws IOException {
        IntegrationTestSpec spec = new IntegrationTestSpec(
                " -mode INDEL" +
                        " -variant " + getToolTestDataDir() + "VQSR.mixedTest.input" +
                        " -vcfOut %s" +
                        " -tranchesFile " + getToolTestDataDir() + "VQSR.mixedTest.tranches" +
                        " -recalFile " + getToolTestDataDir() + "VQSR.mixedTest.recal",
                Arrays.asList(getToolTestDataDir() + "expected.VQSR.mixedTest.noRanges.modeINDEL.vcf")
        );
        spec.executeTest("testApplyRecalibrationIndelNoIntervals", this);
    }

    @Test
    public void testApplyRecalibrationSnpAndIndelTogetherNoIntervals() throws IOException {
        IntegrationTestSpec spec = new IntegrationTestSpec(
                        " -mode BOTH" +
                        " -variant " + getToolTestDataDir() + "VQSR.mixedTest.input" +
                        " -vcfOut %s" +
                        " -tranchesFile " + getToolTestDataDir() + "VQSR.mixedTest.tranches" +
                        " -recalFile " + getToolTestDataDir() + "VQSR.mixedTest.recal",
                Arrays.asList(getToolTestDataDir() + "expected.VQSR.mixedTest.noRanges.vcf")
        );
        spec.executeTest("testApplyRecalibrationSnpAndIndelTogetherNoIntervals", this);
    }

    @Test
    public void testApplyRecalibrationSnpAndIndelTogetherNoIntervalsSomeTrainingSites() throws IOException {
        IntegrationTestSpec spec = new IntegrationTestSpec(
                " -mode BOTH" +
                        " -variant " + getToolTestDataDir() + "VQSR.mixedTest.input" +
                        " -vcfOut %s" +
                        " -tranchesFile " + getToolTestDataDir() + "VQSR.mixedTest.tranches" +
                        " -recalFile " + getToolTestDataDir() + "VQSR.mixedTest.TrainingSites.recal",
                Arrays.asList(getToolTestDataDir() + "expected.VQSR.mixedTest.TrainingSites.noRanges.vcf")
        );
        spec.executeTest("testApplyRecalibrationSnpAndIndelTogetherNoIntervals", this);
    }

    @Test
    public void testApplyRecalibrationSnpAndIndelTogether() throws IOException {
        IntegrationTestSpec spec = new IntegrationTestSpec(
                " -R " + hg19_chr1_1M_Reference +
                        " -L 1:900100-900400" +   //small slice
                        " -mode BOTH" +
                        " -variant " + getToolTestDataDir() + "VQSR.mixedTest.hg19_chr1_1M.input" +
                        " -vcfOut %s" +
                        " -tranchesFile " + getToolTestDataDir() + "VQSR.mixedTest.hg19_chr1_1M.tranches" +
                        " -recalFile " + getToolTestDataDir() + "VQSR.mixedTest.hg19_chr1_1M.recal",
                Arrays.asList(getToolTestDataDir() + "expected.VQSR.mixedTest.hg19_chr1_1M.vcf")
                );
        spec.executeTest("testApplyRecalibrationSnpAndIndelTogether", this);
    }
//
//    @Test(enabled = true)
//    public void testApplyRecalibrationSnpAndIndelTogetherExcludeFiltered() throws Exception {
//        final String base = "-R " + b37KGReference +
//                " -T ApplyRecalibration" +
//                " -L 20:1000100-1000500" +
//                " -mode BOTH" +
//                " --excludeFiltered -ts_filter_level 90.0" +
//                " --no_cmdline_in_header" +
//                " -input " + privateTestDir + "VQSR.mixedTest.input" +
//                " -o %s" +
//                " -tranchesFile " + privateTestDir + "VQSR.mixedTest.tranches" +
//                " -recalFile " + privateTestDir + "VQSR.mixedTest.recal";
//
//        final WalkerTestSpec spec = new WalkerTestSpec(base, 1, Arrays.asList(""));
//        final List<File> outputFiles = executeTest("testApplyRecalibrationSnpAndIndelTogether", spec).getFirst();
//        final File VCF = outputFiles.get(0);
//        for( final VariantContext VC : VCIterable.readAllVCs(VCF, new VCFCodec()).getSecond() ) {
//            if( VC != null ) {
//                Assert.assertTrue(VC.isNotFiltered()); // there should only be unfiltered records in the output VCF file
//            }
//        }
//    }
}
