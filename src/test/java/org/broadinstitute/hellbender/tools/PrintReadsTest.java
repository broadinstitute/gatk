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
import java.nio.file.Files;
import java.util.Arrays;

public class PrintReadsTest extends CommandLineProgramTest{

    private static final File TEST_DATA_DIR = new File("src/test/resources/org/broadinstitute/hellbender/tools/");

    @Override
    public String getCommandLineProgramName() {
        return PrintReads.class.getSimpleName();
    }


    @Test(dataProvider="testingData")
    public void testFileToFile(String fileIn, String extOut) throws Exception {
        String samFile= fileIn;
        final File outFile = File.createTempFile(samFile + ".", extOut);
        outFile.deleteOnExit();
        File ORIG_BAM = new File(TEST_DATA_DIR, samFile);
        final String[] args = new String[]{
                "INPUT=" + ORIG_BAM.getAbsolutePath(),
                "OUTPUT=" + outFile.getAbsolutePath()
        };
        Assert.assertEquals(runCommandLine(args), 0);

        final CompareSAMs compareSAMs = new CompareSAMs();
        compareSAMs.samFiles = Arrays.asList(ORIG_BAM, outFile);
        compareSAMs.doWork();
        Assert.assertTrue(compareSAMs.areEqual());
    }

    @DataProvider(name="testingData")
    public Object[][] testingData() {
        return new String[][]{
                {"print_reads.sam", ".sam"},
                {"print_reads.sam", ".bam"},
                {"print_reads.bam", ".sam"},
                {"print_reads.bam", ".bam"},
        };
    }

}