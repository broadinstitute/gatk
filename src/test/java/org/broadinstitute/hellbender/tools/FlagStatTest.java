/*
 * The MIT License
 *
 * Copyright (c) 2012 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */
package org.broadinstitute.hellbender.tools;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;

public class FlagStatTest extends CommandLineProgramTest{

    @Test(dataProvider = "filenames")
    public void testSamCount(String fileIn) throws Exception {
        final File ORIG_BAM = new File(getTestDataDir(), fileIn);
        final String[] args = new String[]{
                "--input" , ORIG_BAM.getAbsolutePath(),
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
        return new String[][]{
                {"flag_stat.sam"},
                {"flag_stat.bam"},
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
