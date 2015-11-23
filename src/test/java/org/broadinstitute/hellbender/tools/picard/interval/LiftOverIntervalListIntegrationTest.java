package org.broadinstitute.hellbender.tools.picard.interval;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.test.IntegrationTestSpec;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

public final class LiftOverIntervalListIntegrationTest extends CommandLineProgramTest {
    private static final File TEST_DATA_DIR = new File(getTestDataDir(), "picard/interval/" + LiftOverIntervalList.class.getSimpleName());

    @Test
    public void test() throws IOException {
        final File inputInterval = new File(TEST_DATA_DIR, "Broad.human.exome.b37.200lines.interval_list");
        final File inputChain = new File(TEST_DATA_DIR, "b37tohg18.chr1.chain");
        final File inputDict = new File(TEST_DATA_DIR, "Homo_sapiens_assembly18.dict");
        final File expectedFile = new File(TEST_DATA_DIR, "Broad.human.exome.b37.200lines.hg18.interval_list");
        final File outfile = BaseTest.createTempFile("test_Broad.human.exome.b37.200lines.", ".interval_list");
        final String[] args = new String[]{
                "--input", inputInterval.getAbsolutePath(),
                "--CHAIN", inputChain.getAbsolutePath(),
                "--SEQUENCE_DICTIONARY", inputDict.getAbsolutePath(),
                "--output", outfile.getAbsolutePath(),
        };
        final Object code = runCommandLine(args);
        Assert.assertEquals(code, 0);
        IntegrationTestSpec.assertEqualTextFiles(outfile, expectedFile);
    }

    @Test
    public void testBadness() {
        final File inputInterval = new File(TEST_DATA_DIR, "Broad.human.exome.b37.200lines.badness.interval_list");
        final File inputChain = new File(TEST_DATA_DIR, "b37tohg18.chr1.chain");
        final File inputDict = new File(TEST_DATA_DIR, "Homo_sapiens_assembly18.dict");
        final File outfile = BaseTest.createTempFile("test_Broad.human.exome.b37.200lines.", ".interval_list");
        final String[] args = new String[]{
                "--input", inputInterval.getAbsolutePath(),
                "--CHAIN", inputChain.getAbsolutePath(),
                "--SEQUENCE_DICTIONARY", inputDict.getAbsolutePath(),
                "--output", outfile.getAbsolutePath(),
        };
        final Object code = runCommandLine(args);
        Assert.assertEquals(code, 1);
    }

}
