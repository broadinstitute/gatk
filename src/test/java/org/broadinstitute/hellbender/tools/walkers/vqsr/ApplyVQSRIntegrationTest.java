package org.broadinstitute.hellbender.tools.walkers.vqsr;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.tools.IntegrationTestSpec;
import org.testng.annotations.Test;

import java.io.IOException;
import java.util.Arrays;

public class ApplyVQSRIntegrationTest extends CommandLineProgramTest{

    @Test
    public void testApplyRecalibrationSnpAndIndelTogetherNoIntervals() throws IOException {
        IntegrationTestSpec spec = new IntegrationTestSpec(
                        " -mode BOTH" +
                        " -variant " + getToolTestDataDir() + "VQSR.mixedTest.input" +
                        " -vcfOut %s" +
                        " -tranchesFile " + getToolTestDataDir() + "VQSR.mixedTest.tranches" +
                        " -recalFile " + getToolTestDataDir() + "VQSR.mixedTest.recal",
                Arrays.asList("expected.txt")
        );
        spec.executeTest("testApplyRecalibrationSnpAndIndelTogether", this);
    }

    @Test
    public void testApplyRecalibrationSnpAndIndelTogether() throws IOException {
        IntegrationTestSpec spec = new IntegrationTestSpec(
                        " -L 20:1000100-1000500" +
                        " -mode BOTH" +
                        " -variant " + getToolTestDataDir() + "VQSR.mixedTest.input" +
                        " -vcfOut %s" +
                        " -tranchesFile " + getToolTestDataDir() + "VQSR.mixedTest.tranches" +
                        " -recalFile " + getToolTestDataDir() + "VQSR.mixedTest.recal",
                Arrays.asList("expected.txt")
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
