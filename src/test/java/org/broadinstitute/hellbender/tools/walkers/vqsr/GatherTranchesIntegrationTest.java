package org.broadinstitute.hellbender.tools.walkers.vqsr;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;

/**
 * Created by gauthier on 7/18/17.
 */
public class GatherTranchesIntegrationTest extends CommandLineProgramTest {

    private static final String testDir = GATKBaseTest.publicTestDir + "/large/VQSR/";

    @Test
    public void testCombine2Shards() throws Exception {
        final File recal1 = new File(testDir + "snpTranches.scattered.txt");  //this is the output of VariantRecalibratorIntegrationTest.testVariantRecalibratorSNPscattered
        final File recal2 = new File(testDir + "snpTranches.scattered.2.txt"); //this is a copy of the above

        final File recal_original = new File(testDir + "expected/snpTranches.gathered.txt");

        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addRaw("--input");
        args.addRaw(recal1.getAbsolutePath());
        args.addRaw("--input");
        args.addRaw(recal2.getAbsolutePath());
        args.add("mode", "SNP");

        final File outFile = GATKBaseTest.createTempFile("gatheredTranches", ".txt");
        args.addOutput(outFile);
        final Object res = this.runCommandLine(args.getArgsArray());
        Assert.assertEquals(res, 0);
        IntegrationTestSpec.assertEqualTextFiles(outFile, recal_original);
    }

    @Test
    public void testCombine2IndelTranches() throws Exception {
        final File tranches1 = new File(testDir + "indels.0.tranches");
        final File tranches2 = new File(testDir + "indels.1.tranches");

        final File recal_original = new File(testDir + "expected/indels.gathered.tranches");

        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addRaw("--input");
        args.addRaw(tranches1.getAbsolutePath());
        args.addRaw("--input");
        args.addRaw(tranches2.getAbsolutePath());
        args.add("mode", "INDEL");

        final File outFile = GATKBaseTest.createTempFile("gatheredTranches", ".txt");
        args.addOutput(outFile);
        final Object res = this.runCommandLine(args.getArgsArray());
        Assert.assertEquals(res, 0);
        IntegrationTestSpec.assertEqualTextFiles(outFile, recal_original);
    }
}