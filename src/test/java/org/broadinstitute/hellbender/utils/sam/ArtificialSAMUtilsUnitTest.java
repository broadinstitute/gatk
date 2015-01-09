/*
* Copyright (c) 2012 The Broad Institute
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

package org.broadinstitute.hellbender.utils.sam;

import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.utils.BaseTest;
import org.broadinstitute.hellbender.utils.iterators.GATKSAMIterator;
import org.testng.annotations.Test;

import static org.testng.Assert.*;

public class ArtificialSAMUtilsUnitTest extends BaseTest {


    @Test
    public void basicReadIteratorTest() {
        GATKSAMIterator iter = ArtificialSAMUtils.mappedReadIterator(1, 100, 100);
        int count = 0;
        while (iter.hasNext()) {
            SAMRecord rec = iter.next();
            count++;
        }
        assertEquals(count, 100 * 100);
    }

    @Test
    public void tenPerChromosome() {
        GATKSAMIterator iter = ArtificialSAMUtils.mappedReadIterator(1, 100, 10);
        int count = 0;
        while (iter.hasNext()) {
            SAMRecord rec = iter.next();

            assertEquals(Integer.valueOf(Math.round(count / 10)), rec.getReferenceIndex());
            count++;
        }
        assertEquals(count, 100 * 10);
    }

    @Test
    public void onePerChromosome() {
        GATKSAMIterator iter = ArtificialSAMUtils.mappedReadIterator(1, 100, 1);
        int count = 0;
        while (iter.hasNext()) {
            SAMRecord rec = iter.next();

            assertEquals(Integer.valueOf(count), rec.getReferenceIndex());
            count++;
        }
        assertEquals(count, 100 * 1);
    }

    @Test
    public void basicUnmappedIteratorTest() {
        GATKSAMIterator iter = ArtificialSAMUtils.mappedAndUnmappedReadIterator(1, 100, 100, 1000);
        int count = 0;
        for (int x = 0; x < (100* 100); x++ ) {
            if (!iter.hasNext()) {
                fail ("we didn't get the expected number of reads");
            }
            SAMRecord rec = iter.next();
            assertTrue(rec.getReferenceIndex() >= 0);
            count++;
        }
        assertEquals(100 * 100, count);

        // now we should have 1000 unmapped reads
        count = 0;
        while (iter.hasNext()) {
            SAMRecord rec = iter.next();
            assertTrue(rec.getReferenceIndex() < 0);
            count++;
        }
        assertEquals(count, 1000);
    }


}
