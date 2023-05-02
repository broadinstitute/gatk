package org.broadinstitute.hellbender.tools.walkers.gvs;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;

public class GvsPgenExtractorIntegrationTest extends CommandLineProgramTest {

    @DataProvider(name = "argumentSets")
    public Object[][] getArgumentSets() {
        return new Object[][]{
            {
                new ArgumentsBuilder()
                        .addVCF(getToolTestDataDir() + "CEUtrioTest.vcf")
            }
        };
    }

    @Test(dataProvider = "argumentSets")
    public void testGvsPgenExtractor(
        ArgumentsBuilder commandLineArgs
    ) {
        // Run each set of command line arguments to create a tmp pgen
        final File outputDir = createTempDir("testGvsPgenExtractor");
        final File outputFile = new File(outputDir + "/extracted.pgen");
        final ArgumentsBuilder args = commandLineArgs.addOutput(outputFile);
        runCommandLine(args, GvsPgenExtractor.class.getSimpleName());
    }

}
