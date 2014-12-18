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

import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;

public class FlagStatTest {
    @Test
    public void testSamCount() throws Exception {
        final FlagStat stats = new FlagStat();
        stats.INPUT = new File("src/test/resources/org/broadinstitute/hellbender/tools/flag_stat.sam");

        final FlagStat.FlagStatus l = stats.countReads();
        Assert.assertEquals(l.readCount, 19);
        Assert.assertEquals(l.QC_failure, 2);
        Assert.assertEquals(l.duplicates, 0);
        Assert.assertEquals(l.mapped, 11);
        Assert.assertEquals(l.paired_in_sequencing, 19);
        Assert.assertEquals(l.read1, 9);
        Assert.assertEquals(l.read2, 10);
        Assert.assertEquals(l.properly_paired, 5);
        Assert.assertEquals(l.with_itself_and_mate_mapped, 5);
        Assert.assertEquals(l.singletons, 6);
        Assert.assertEquals(l.with_mate_mapped_to_a_different_chr, 0);
        Assert.assertEquals(l.with_mate_mapped_to_a_different_chr_maq_greaterequal_than_5, 0);
    }

    @Test
    public void testBamCount() throws Exception {
        final FlagStat stats = new FlagStat();
        stats.INPUT = new File("src/test/resources/org/broadinstitute/hellbender/tools/flag_stat.bam");

        final FlagStat.FlagStatus l = stats.countReads();
        Assert.assertEquals(l.readCount, 19);
        Assert.assertEquals(l.QC_failure, 2);
        Assert.assertEquals(l.duplicates, 0);
        Assert.assertEquals(l.mapped, 11);
        Assert.assertEquals(l.paired_in_sequencing, 19);
        Assert.assertEquals(l.read1, 9);
        Assert.assertEquals(l.read2, 10);
        Assert.assertEquals(l.properly_paired, 5);
        Assert.assertEquals(l.with_itself_and_mate_mapped, 5);
        Assert.assertEquals(l.singletons, 6);
        Assert.assertEquals(l.with_mate_mapped_to_a_different_chr, 0);
        Assert.assertEquals(l.with_mate_mapped_to_a_different_chr_maq_greaterequal_than_5, 0);
    }
}
