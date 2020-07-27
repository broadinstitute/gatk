package org.broadinstitute.hellbender.tools.walkers;

import static org.testng.Assert.*;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.*;


public class MatchWithReferencePanelTest extends CommandLineProgramTest {

    @Override
    public String getTestedClassName() {
        return MatchWithReferencePanel.class.getSimpleName();
    }


    @Test(dataProvider = "matchWithReferencePanelInputs")
    public void testMatchWithReferencePanel(final File fileIn, final File panelFile) throws Exception {


        final File outputVCF = createTempFile("testMatchWithReferencePanel", ".vcf");
        final File expectedVCF = getTestFile("matchWithReferencePanel.expected.vcf");

        final ArgumentsBuilder ab = new ArgumentsBuilder();
        ab.addVCF(fileIn);
        ab.add("panel", panelFile);
       // Arrays.asList("-panel",
        ab.addOutput(outputVCF);

        runCommandLine(ab.getArgsArray());

        // Test for an exact match against past results
        IntegrationTestSpec.assertEqualTextFiles(outputVCF, expectedVCF);

    }

    @DataProvider(name="matchWithReferencePanelInputs")
    public Object[][] generateInputs() {
        return new Object[][]{
                {getTestFile("testArray.vcf"), getTestFile("testPanel.vcf")},

        };
    }
}
