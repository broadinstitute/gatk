package org.broadinstitute.hellbender.tools.walkers.vqsr;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;

/**
 * Created by gauthier on 7/18/17.
 */
public class GatherTranchesIntegrationTest extends CommandLineProgramTest {

    private static final String testDir = GATKBaseTest.publicTestDir + "/large/VQSR/";

    @DataProvider(name = "testInputs")
    public Object[][] getTestInputs () {
        return new Object[][]{
                {Arrays.asList(new File(testDir + "snpTranches.scattered.txt"), new File(testDir + "snpTranches.scattered.txt")),
                        new File(testDir + "expected/snpTranches.gathered.txt"), "SNP"},

                {Arrays.asList(new File(testDir + "indels.0.tranches"), new File(testDir + "indels.1.tranches")),
                        new File(testDir + "expected/indels.gathered.tranches"), "INDEL"},

                {Arrays.asList(new File(testDir + "test-single-giant-input-snps.tranches")),
                        new File(testDir + "expected/singleOverflow.tranches"), "SNP"},

                {Arrays.asList(new File(testDir + "test-very-large-one-snps.tranches"), new File(testDir + "test-very-large-two-snps.tranches")),
                        new File(testDir + "expected/testSummedOverflow.tranches"), "SNP"}
        };
    }

    @Test (dataProvider = "testInputs")
    public void testGatherTranches(List<File> inputs, File expected, String mode) throws IOException {
        final ArgumentsBuilder args = new ArgumentsBuilder();
        for (File inFile : inputs) {
            args.addRaw("--input");
            args.addRaw(inFile);
        }
        args.add("mode", mode);

        final File outFile = GATKBaseTest.createTempFile("gatheredTranches", ".txt");
        args.addOutput(outFile);
        final Object res = this.runCommandLine(args.getArgsArray());
        Assert.assertEquals(res, 0);
        IntegrationTestSpec.assertEqualTextFiles(outFile, expected);
    }
}