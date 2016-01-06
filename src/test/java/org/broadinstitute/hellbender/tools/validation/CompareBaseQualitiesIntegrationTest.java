package org.broadinstitute.hellbender.tools.validation;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.test.IntegrationTestSpec;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;

public class CompareBaseQualitiesIntegrationTest extends CommandLineProgramTest {
    @Test
    public void identicalBamTest() throws Exception {
        // This test verifies that if two bams are identical that they produce zero diffs.
        final String resourceDir = getTestDataDir() + "/validation/";

        final File firstBam = new File(resourceDir, "CEUTrio.HiSeq.WGS.b37.ch20.1m-1m1k.NA12878.bam");
        final File secondBam = new File(resourceDir, "CEUTrio.HiSeq.WGS.b37.ch20.1m-1m1k.NA12878.bam");

        ArgumentsBuilder args = new ArgumentsBuilder();
        args.add(firstBam.getCanonicalPath());
        args.add(secondBam.getCanonicalPath());
        args.add("--throwOnDiff true");
        args.add("--VALIDATION_STRINGENCY SILENT");
        this.runCommandLine(args.getArgsArray());
    }

    @DataProvider(name = "CompareBasesProvider")
    public Object[][] makeCompareBasesProvider() {
        final String resourceDir = getTestDataDir() + "/validation/";
        final File outFile = BaseTest.createTempFile(getTestDataDir() + "temp.diff", "txt");
        final File firstBam = new File(resourceDir, "single.read.bam");
        final File secondBam = new File(resourceDir, "another.single.read.bam");
        return new Object[][]{
                {firstBam, firstBam, outFile, "single.read.qual.diff.txt"},
                {firstBam, secondBam, outFile, "two.reads.qual.diff.txt"},
        };

    }

    @Test(dataProvider = "CompareBasesProvider")
    public void singleReadDiffTest(File firstBam, File secondBam, File outFile, String diffFile) throws Exception {
        final String resourceDir = getTestDataDir() + "/validation/";

        ArgumentsBuilder args = new ArgumentsBuilder();
        args.add(firstBam.getCanonicalPath());
        args.add(secondBam.getCanonicalPath());
        args.add("--" + "O");
        args.add(outFile.getCanonicalPath());

        this.runCommandLine(args.getArgsArray());


        final File expected = new File(resourceDir, diffFile);
        IntegrationTestSpec.assertEqualTextFiles(outFile, expected);
    }
}