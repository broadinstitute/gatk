/*
* Copyright 2012-2015 Broad Institute, Inc.
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.gatk.utils.sam;

import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMRecord;
import org.broadinstitute.gatk.utils.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

/**
 * Created with IntelliJ IDEA.
 * User: depristo
 * Date: 1/15/13
 * Time: 3:49 PM
 * To change this template use File | Settings | File Templates.
 */
public class ArtificialBAMBuilderUnitTest extends BaseTest {
    @DataProvider(name = "ArtificialBAMBuilderUnitTestProvider")
    public Object[][] makeArtificialBAMBuilderUnitTestProvider() {
        final List<Object[]> tests = new LinkedList<Object[]>();

        final List<Integer> starts = Arrays.asList(
                1, // very start of the chromosome
                ArtificialBAMBuilder.BAM_SHARD_SIZE - 100, // right before the shard boundary
                ArtificialBAMBuilder.BAM_SHARD_SIZE + 100 // right after the shard boundary
        );

        for ( final int readLength : Arrays.asList(10, 20) ) {
            for ( final int skips : Arrays.asList(0, 1, 10) ) {
                for ( final int start : starts ) {
                    for ( final int nSamples : Arrays.asList(1, 2) ) {
                        for ( final int nReadsPerLocus : Arrays.asList(1, 10) ) {
                            for ( final int nLoci : Arrays.asList(10, 100, 1000) ) {
                                final ArtificialBAMBuilder bamBuilder = new ArtificialBAMBuilder(nReadsPerLocus, nLoci);
                                bamBuilder.setReadLength(readLength);
                                bamBuilder.setSkipNLoci(skips);
                                bamBuilder.setAlignmentStart(start);
                                bamBuilder.createAndSetHeader(nSamples);
                                tests.add(new Object[]{bamBuilder, readLength, skips, start, nSamples, nReadsPerLocus, nLoci});
                            }
                        }
                    }
                }
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "ArtificialBAMBuilderUnitTestProvider")
    public void testBamProvider(final ArtificialBAMBuilder bamBuilder, int readLength, int skips, int start, int nSamples, int nReadsPerLocus, int nLoci) {
        Assert.assertEquals(bamBuilder.getReadLength(), readLength);
        Assert.assertEquals(bamBuilder.getSkipNLoci(), skips);
        Assert.assertEquals(bamBuilder.getAlignmentStart(), start);
        Assert.assertEquals(bamBuilder.getNSamples(), nSamples);
        Assert.assertEquals(bamBuilder.getnReadsPerLocus(), nReadsPerLocus);
        Assert.assertEquals(bamBuilder.getnLoci(), nLoci);

        final List<GATKSAMRecord> reads = bamBuilder.makeReads();
        Assert.assertEquals(reads.size(), bamBuilder.expectedNumberOfReads());
        for ( final GATKSAMRecord read : reads ) {
            assertGoodRead(read, bamBuilder);
        }

        final File bam = bamBuilder.makeTemporarilyBAMFile();
        final SAMFileReader reader = new SAMFileReader(bam);
        Assert.assertTrue(reader.hasIndex());
        final Iterator<SAMRecord> bamIt = reader.iterator();
        int nReadsFromBam = 0;
        int lastStart = -1;
        while ( bamIt.hasNext() ) {
            final SAMRecord read = bamIt.next();
            assertGoodRead(read, bamBuilder);
            nReadsFromBam++;
            Assert.assertTrue(read.getAlignmentStart() >= lastStart);
            lastStart = read.getAlignmentStart();
        }
        Assert.assertEquals(nReadsFromBam, bamBuilder.expectedNumberOfReads());
    }

    private void assertGoodRead(final SAMRecord read, final ArtificialBAMBuilder bamBuilder) {
        Assert.assertEquals(read.getReadLength(), bamBuilder.getReadLength());
        Assert.assertEquals(read.getReadBases().length, bamBuilder.getReadLength());
        Assert.assertEquals(read.getBaseQualities().length, bamBuilder.getReadLength());
        Assert.assertTrue(read.getAlignmentStart() >= bamBuilder.getAlignmentStart());
        Assert.assertNotNull(read.getReadGroup());
    }
}


