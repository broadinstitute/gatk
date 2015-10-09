package org.broadinstitute.hellbender.tools;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.test.SamAssertionUtils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
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

    @Test(dataProvider = "clipOptions")
    public void testClipper(String inBam, String option, String optAbrv) throws IOException {
        final String tmpName = BaseTest.createTempFile(inBam + "." + optAbrv, ".tmp").getAbsolutePath();
        final String outName = BaseTest.createTempFile(inBam + "." + optAbrv, ".bam").getAbsolutePath();

        final File outFileBam = new File(outName);
        final File outFileTmp = new File(tmpName);
        final List<String> args = new ArrayList<>();
        args.addAll(Arrays.<String>asList(
                "--input", new File(getTestDataDir(), inBam + ".bam").getAbsolutePath(),
                "-os", tmpName,
                "--output", outName
        ));
        args.addAll(Arrays.asList(option.split("\\s+")));

        final File expectedOutBam = new File(getTestDataDir(), "expected." + inBam + "." + optAbrv + ".bam");
        final File expectedTmp = new File(getTestDataDir(), "expected." + inBam + "." + optAbrv + ".tmp");
        final ClipReads.ClippingData res = (ClipReads.ClippingData)this.runCommandLine(args);
        System.out.println(res);

        Assert.assertTrue(expectedOutBam.exists(), "expected output read file exists " + expectedOutBam.getAbsolutePath());
        Assert.assertTrue(outFileBam.exists(), "actual output read file exists " + outFileBam.getAbsolutePath());

        Assert.assertTrue(expectedTmp.exists(), "expected output stat file exists " + expectedTmp.getAbsolutePath());
        Assert.assertTrue(outFileTmp.exists(), "actual output stat file exists " + outFileTmp.getAbsolutePath());
        SamAssertionUtils.assertSamsEqual(expectedOutBam, outFileBam);

        List<String> actualLines = new XReadLines(new File(tmpName)).readLines();
        List<String> expectedLines = new XReadLines(expectedTmp).readLines();
        Assert.assertEquals(actualLines.toString(), expectedLines.toString());
    }

    @DataProvider(name="clipOptions")
    public Object[][] clipOptions() {
        final String b1 = "clippingReadsTest.withRG.hg19";
        return new String[][]{
                {b1, "-QT 0", "QT_0"},
                {b1, "-QT 2", "QT_2"},
                {b1, "-QT 10", "QT_10"},
                {b1, "-QT 20", "QT_20"},
                {b1, "-CT 1-5", "CT_15"},
                {b1, "-CT 1-5,11-15", "CT_15_1115"},
                {b1, "-X CCCCC", "X_CCCCC"},
                {b1, "-QT 10 -CR WRITE_NS", "QT_10_CR_WRITE_NS"},
                {b1, "-QT 10 -CR WRITE_Q0S", "QT_10_CR_WRITE_Q0S"},
                {b1, "-QT 10 -CR SOFTCLIP_BASES","QT_10_CR_SOFTCLIP_BASES"} ,
                {b1, "-XF " + BaseTest.publicTestDir + "seqsToClip.fasta", "XF"},
                {b1, "-QT 10 -CT 1-5 -X CCCCC -XF " + BaseTest.publicTestDir + "seqsToClip.fasta", "QT_10_CT_15_X_CCCCC_XF"},
        };
    }
}
