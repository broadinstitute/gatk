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
        final File recal2 = new File(testDir + "snpTranches.scattered.2.txt"); //this is the output of running the equivalent commandline as above outside of the integration test :-/

        final File recal_original = new File(testDir + "expected/snpTranches.gathered.txt");

        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("--input");
        args.add(recal1.getAbsolutePath());
        args.add("--input");
        args.add(recal2.getAbsolutePath());

        final File outFile = GATKBaseTest.createTempFile("gatheredTranches", ".txt");
        args.addOutput(outFile);
        final Object res = this.runCommandLine(args.getArgsArray());
        Assert.assertEquals(res, 0);
        IntegrationTestSpec.assertEqualTextFiles(outFile, recal_original);
    }


}