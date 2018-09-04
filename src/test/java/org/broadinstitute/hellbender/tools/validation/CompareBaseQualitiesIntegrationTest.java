package org.broadinstitute.hellbender.tools.validation;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Objects;
import java.util.Arrays;
import java.util.List;

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
        args.add("--throw-on-diff true");
        args.add("--VALIDATION_STRINGENCY SILENT");
        final Object result = this.runCommandLine(args);
        Assert.assertEquals(result, 0);
    }

    @DataProvider(name = "CompareBasesProvider")
    public Object[][] makeCompareBasesProvider() {
        final String resourceDir = getTestDataDir() + "/validation/";
        final File outFile = GATKBaseTest.createTempFile(getTestDataDir() + "temp.diff", "txt");
        final File firstBam = new File(resourceDir, "single.read.bam");
        final File secondBam = new File(resourceDir, "another.single.read.bam");
        final File firstCram = new File(resourceDir, "single.read.cram");
        final File secondCram = new File(resourceDir, "another.single.read.cram");
        final File referenceFile = new File(b37_reference_20_21);

        final List<Integer> sq = Arrays.asList(10, 20, 30, 40);
        return new Object[][]{
                {firstBam, firstBam, null, outFile, sq, "single.read.qual.diff.txt"},
                {firstBam, secondBam, null, outFile, sq, "two.reads.qual.diff.txt"},
                {firstCram, secondCram, referenceFile, outFile, sq, "two.reads.qual.diff.txt"},
                {firstBam, secondCram, referenceFile, outFile, sq, "two.reads.qual.diff.txt"},
                {firstCram, secondBam, referenceFile, outFile, sq, "two.reads.qual.diff.txt"},
        };

    }

    @Test(dataProvider = "CompareBasesProvider")
    public void singleReadDiffTest(File firstBam, File secondBam, File referenceFile, File outFile, List<Integer> staticQuantizationQuals, String diffFile) throws Exception {
        final String resourceDir = getTestDataDir() + "/validation/";

        ArgumentsBuilder args = new ArgumentsBuilder();
        args.add(firstBam.getCanonicalPath());
        args.add(secondBam.getCanonicalPath());
        args.add("--" + "O");
        args.add(outFile.getCanonicalPath());
        if (null != referenceFile) {
            args.add("-R");
            args.add(referenceFile.getAbsolutePath());
        }
        if (staticQuantizationQuals != null && !staticQuantizationQuals.isEmpty()){
            for (int sq : staticQuantizationQuals){
                args.add("--" + CompareBaseQualities.STATIC_QUANTIZED_QUALS_LONG_NAME);
                args.add(sq);
            }
        }

        final Object result = this.runCommandLine(args);
        if (Objects.equals(firstBam, secondBam)) {
            Assert.assertEquals(result, 0);
        } else {
            Assert.assertNotEquals(result, 0);
        }

        final File expected = new File(resourceDir, diffFile);
        IntegrationTestSpec.assertEqualTextFiles(outFile, expected);
    }
}