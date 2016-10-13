package org.broadinstitute.hellbender.tools;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;

public final class FlagStatIntegrationTest extends CommandLineProgramTest{

    @Test(dataProvider = "filenames")
    public void testSamCount(String fileIn, boolean useSpark) throws Exception {
        final File ORIG_BAM = new File(getTestDataDir(), fileIn);
        final String[] args = new String[]{
                "--input" , ORIG_BAM.getAbsolutePath(),
                "--spark" , Boolean.toString(useSpark),
        };
        final Object res = this.runCommandLine(args);

        final FlagStat.FlagStatus l = new FlagStat.FlagStatus();
        l.readCount = 19;
        l.QC_failure = 2;
        l.duplicates = 0;
        l.mapped = 11;
        l.paired_in_sequencing = 19;
        l.read1 = 9;
        l.read2 = 10;
        l.properly_paired = 5;
        l.with_itself_and_mate_mapped = 5;
        l.singletons = 6;
        l.with_mate_mapped_to_a_different_chr = 0;
        l.with_mate_mapped_to_a_different_chr_maq_greaterequal_than_5 = 0;
        Assert.assertEquals(res, l);
    }


    @DataProvider(name="filenames")
    public Object[][] filenames() {
        return new Object[][]{
                {"flag_stat.sam", false},
                {"flag_stat.sam", true},
                {"flag_stat.bam", false},
                {"flag_stat.bam", true},
        };
    }

    @Test
    public void testEqualFS(){
        FlagStat.FlagStatus l1 = makeFlagStatus();
        FlagStat.FlagStatus l2 = makeFlagStatus();
        Assert.assertEquals(l1, l2);
        Assert.assertEquals(l1.hashCode(), l2.hashCode());
        Assert.assertNotSame(l1, l2);
    }

    @Test
    public void testNonEqualFS(){
        FlagStat.FlagStatus l1 = makeFlagStatus();
        FlagStat.FlagStatus l2 = makeFlagStatus();
        l2.duplicates++;
        Assert.assertNotEquals(l1, l2);
        Assert.assertNotEquals(l1.hashCode(), l2.hashCode());
        Assert.assertNotSame(l1, l2);
    }

    private FlagStat.FlagStatus makeFlagStatus() {
        FlagStat.FlagStatus l = new FlagStat.FlagStatus();
        l.readCount = 19;
        l.QC_failure = 2;
        l.duplicates = 0;
        l.mapped = 11;
        l.paired_in_sequencing = 19;
        l.read1 = 9;
        l.read2 = 10;
        l.properly_paired = 5;
        l.with_itself_and_mate_mapped = 5;
        l.singletons = 6;
        l.with_mate_mapped_to_a_different_chr = 0;
        l.with_mate_mapped_to_a_different_chr_maq_greaterequal_than_5 = 0;
        return l;
    }
}
