package org.broadinstitute.hellbender.tools.walkers.sv;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;

public class EvaluateConcordanceOnNonGenotypedSVVCFIntegrationTest extends CommandLineProgramTest {

    @Override
    public String getTestedClassName() {
        return EvaluateConcordanceOnNonGenotypedSVVCF.class.getSimpleName();
    }

    @Test(groups = "sv")
    public void testRunnable() {

        final ArgumentsBuilder argumentsBuilder = new ArgumentsBuilder();

        // near minimal set of arguments
        argumentsBuilder.addArgument(EvaluateConcordanceOnNonGenotypedSVVCF.TRUTH_VARIANTS_SHORT_NAME, largeFileTestDir + "sv/ConcordanceTest_Truth.vcf");
        argumentsBuilder.addArgument(EvaluateConcordanceOnNonGenotypedSVVCF.NON_CANONICAL_CHROMOSOMES_FILE_SHORT_NAME, toolsTestDir + "spark/sv/integration/inputs/Homo_sapiens_assembly38.kill.alts");

        argumentsBuilder.addVCF( new File(largeFileTestDir + "sv/ConcordanceTest_Eval.vcf") );

        // output to temp location
        final File bedOutPut = GATKBaseTest.createTempFile("sv_concordance", ".bed");
        argumentsBuilder.addOutput(bedOutPut);

        runCommandLine(argumentsBuilder);

        assertionsOnOutput(bedOutPut);
    }

    private void assertionsOnOutput(final File bedFile) {

        // simply asserting that the BED file exists and is not empty, for now
        Assert.assertTrue(bedFile.isFile());
        Assert.assertTrue(bedFile.length() != 0);
    }
}
