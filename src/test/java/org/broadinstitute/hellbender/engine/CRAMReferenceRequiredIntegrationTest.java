package org.broadinstitute.hellbender.engine;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.IntegrationTestSpec;
import org.broadinstitute.hellbender.tools.PrintReads;
import org.broadinstitute.hellbender.utils.read.SamAssertionUtils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Collections;

/**
 * Created by edwardk on 7/6/15.
 */
public class CRAMReferenceRequiredIntegrationTest extends CommandLineProgramTest {

    public static final String SEQDICTVAL_TEST_DIRECTORY = "src/test/resources/org/broadinstitute/hellbender/utils/CRAMReference/";

    @Override
    public String getTestedClassName() {
        return PrintReads.class.getSimpleName();
    }

    @Test
    public void testReferenceGiven() throws IOException {
        String samFile= "refrequired.cram";
        final File outFile = BaseTest.createTempFile(samFile + ".", ".bam");
        File ORIG_BAM = new File(SEQDICTVAL_TEST_DIRECTORY, samFile);
        final String[] args = new String[]{
                "--input" , ORIG_BAM.getPath(),
                "--output", outFile.getPath()
        };
        Assert.assertEquals(runCommandLine(args), null);
        SamAssertionUtils.assertSamsEqual(ORIG_BAM, outFile);
    }
    /*
    public void testReferenceGiven() throws IOException {
        IntegrationTestSpec testSpec = new IntegrationTestSpec(
                " -R " + SEQDICTVAL_TEST_DIRECTORY + "refrequired.fasta" +
                " -I " + SEQDICTVAL_TEST_DIRECTORY + "refrequired.cram" +
                " -O " + SEQDICTVAL_TEST_DIRECTORY + "refreqout.txt",
                Collections.emptyList());
        testSpec.executeTest("testReferenceGiven", this);
    }*/

    @Test
    public void testNoReferenceGiven() throws IOException {
        IntegrationTestSpec testSpec = new IntegrationTestSpec(
                        " -I " + SEQDICTVAL_TEST_DIRECTORY + "refrequired.cram" +
                        " -V " + "TestFeatures:" + SEQDICTVAL_TEST_DIRECTORY + "test.vcf",
                0,
                UserException.class);
        testSpec.executeTest("testNoReferenceGiven", this);
    }
}
