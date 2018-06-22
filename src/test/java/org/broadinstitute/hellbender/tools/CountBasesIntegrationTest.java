package org.broadinstitute.hellbender.tools;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;

public final class CountBasesIntegrationTest extends CommandLineProgramTest {

    @Test(dataProvider = "filenames")
    public void testCountBases(final String fileIn, final String referenceName) throws Exception {
        final File ORIG_BAM = new File(getTestDataDir(), fileIn);
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("--input");
        args.add(ORIG_BAM.getAbsolutePath());
        if (null != referenceName) {
            final File REF = new File(getTestDataDir(), referenceName);
            args.add("-R");
            args.add(REF.getAbsolutePath());
        }

        final Object res = this.runCommandLine(args.getArgsArray());
        Assert.assertEquals(res, 808l);
    }

    @DataProvider(name="filenames")
    public Object[][] filenames() {
        return new String[][]{
                {"count_bases.sam", null},
                {"count_bases.bam", null},
                {"count_bases.cram", "count_bases.fasta"}
        };
    }

}
