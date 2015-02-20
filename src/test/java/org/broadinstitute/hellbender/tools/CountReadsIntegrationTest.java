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
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;

public class CountReadsIntegrationTest extends CommandLineProgramTest {
    @Test(dataProvider = "filenames")
    public void testCountBases(String fileIn) throws Exception {
        final File ORIG_BAM = new File(getTestDataDir(), fileIn);
        final String[] args = new String[]{
                "--input",  ORIG_BAM.getAbsolutePath(),
        };
        final Object res = this.runCommandLine(args);
        Assert.assertEquals(res, 8l);
    }

    @DataProvider(name="filenames")
    public Object[][] filenames() {
        return new String[][]{
                {"count_reads.sam"},
                {"count_reads.bam"},
                {"count_reads_sorted.sam"},
                {"count_reads_sorted.bam"}
        };
    }

    @DataProvider(name="intervals")
    public Object[][] intervals(){
        return new Object[][]{
                new Object[]{"-L chr7:1", 3l},
                new Object[]{"", 8l}, //no intervals specified see all reads including unaligned
                new Object[]{"-L chr7:1-20", 4l},
                new Object[]{"-L chr1", 0l},
                new Object[]{"-L chr1 -L chr7", 7l},
                new Object[]{"-XL chr7", 0l},
                new Object[]{"-XL chr7:2-404", 3l},
                new Object[]{"-L chr7:1-30 -L chr7:10-15 --interval_set_rule INTERSECTION", 3l},
                new Object[]{"-L chr7:1 --interval_padding 19", 4l },
                new Object[]{"-L " + getTestDataDir() + "/chr7_1_20.interval_list", 4l},
                new Object[]{"-L chr7:1-100 -XL chr7:2-100", 3l},
                new Object[]{"-L chr7:1-10 -L chr7:5-10 --interval_padding 10 --interval_set_rule INTERSECTION --XL chr7:21-200", 4l }
        };
    }


    @Test(dataProvider = "intervals")
    public void testCountBasesWithIntervals(String interval_args, long count) throws Exception {
        final File ORIG_BAM = new File(getTestDataDir(), "count_reads_sorted.bam");
        ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("--input");
        args.add(ORIG_BAM.getAbsolutePath());
        args.add(interval_args);

        final Object res = this.runCommandLine(args.getArgsArray());
        Assert.assertEquals(res, count);
    }


}
