package org.broadinstitute.hellbender.tools;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.test.SamAssertionUtils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.test.TestResources;
import org.broadinstitute.hellbender.utils.text.XReadLines;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public final class ClipReadsIntegrationTest extends CommandLineProgramTest {

    private File localTestData = new File(getTestDataDir(), "ClipReads");

    @Test(dataProvider = "clipOptions")
    public void testClipper(String inBam, String reference, String extension, String option, String optAbrv, boolean doStats) throws IOException {
        final String tmpBAMOutName = BaseTest.createTempFile(inBam + "." + optAbrv, extension).getAbsolutePath();
        String tmpStatOutName = null;
        if (doStats) {
            tmpStatOutName= BaseTest.createTempFile(inBam + "." + optAbrv, ".tmp").getAbsolutePath();
        }

        final List<String> args = new ArrayList<>();
        args.addAll(Arrays.<String>asList(
                "--input", new File(localTestData, inBam + extension).getAbsolutePath(),
                "--output", tmpBAMOutName
        ));
        if (doStats) {
            args.addAll(Arrays.<String>asList("-os", tmpStatOutName));
        }

        File referenceFile = null;
        if (reference != null) {
            referenceFile = new File(getTestDataDir(),reference);
            args.add("-R");
            args.add(referenceFile.getAbsolutePath());
        }
        args.addAll(Arrays.asList(option.split("\\s+")));

        final ClipReads.ClippingData res = (ClipReads.ClippingData)this.runCommandLine(args);
        System.out.println(res);

        final File outFileBam = new File(tmpBAMOutName);
        final File expectedOutBam = new File(localTestData, "expected." + inBam + "." + optAbrv + extension);
        Assert.assertTrue(expectedOutBam.exists(), "expected output read file exists " + expectedOutBam.getAbsolutePath());
        Assert.assertTrue(outFileBam.exists(), "actual output read file exists " + outFileBam.getAbsolutePath());
        SamAssertionUtils.assertSamsEqual(expectedOutBam, outFileBam, referenceFile);

        if (doStats) {
            final File outFileStat = new File(tmpStatOutName);
            final File expectedOutStat = new File(localTestData, "expected." + inBam + "." + optAbrv + ".tmp");
            Assert.assertTrue(expectedOutStat.exists(), "expected output stat file exists " + expectedOutStat.getAbsolutePath());
            Assert.assertTrue(outFileStat.exists(), "actual output stat file exists " + outFileStat.getAbsolutePath());
            List<String> actualLines = new XReadLines(new File(tmpStatOutName)).readLines();
            List<String> expectedLines = new XReadLines(expectedOutStat).readLines();
            Assert.assertEquals(actualLines.toString(), expectedLines.toString());
        }
    }

    @DataProvider(name="clipOptions")
    public Object[][] clipOptions() {
        final String b1 = "clippingReadsTest.withRG.hg19";
        final String cramFile = "clippingReadsTestCRAM";
        final String referenceFile = "valid.fasta"; // lives in tools test folder
        return new Object[][]{
                {b1, null, ".bam", "-QT 0", "QT_0", true},
                {b1, null, ".bam", "-QT 0", "QT_0", false},
                {b1, null, ".bam", "-QT 2", "QT_2", true},
                {b1, null, ".bam", "-QT 10", "QT_10", true},
                {b1, null, ".bam", "-QT 20", "QT_20", true},
                {b1, null, ".bam", "-CT 1-5", "CT_15", true},
                {b1, null, ".bam", "-CT 1-5,11-15", "CT_15_1115", true},
                {b1, null, ".bam", "-X CCCCC", "X_CCCCC", true},
                {b1, null, ".bam", "-X CCCCC", "X_CCCCC", false},
                {b1, null, ".bam", "-QT 10 -CR WRITE_NS", "QT_10_CR_WRITE_NS", true},
                {b1, null, ".bam", "-QT 10 -CR WRITE_Q0S", "QT_10_CR_WRITE_Q0S", true},
                {b1, null, ".bam", "-QT 10 -CR SOFTCLIP_BASES","QT_10_CR_SOFTCLIP_BASES", true} ,
                {b1, null, ".bam", "-XF " + TestResources.publicTestDir + "seqsToClip.fasta", "XF", true},
                {b1, null, ".bam", "-QT 10 -CT 1-5 -X CCCCC -XF " + TestResources.publicTestDir + "seqsToClip.fasta", "QT_10_CT_15_X_CCCCC_XF", true},
                {cramFile, referenceFile, ".cram", "-QT 10", "QT_10", true},
                {cramFile, referenceFile, ".cram", "-QT 10", "QT_10", false},
        };
    }
}
