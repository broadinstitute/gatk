package org.broadinstitute.hellbender.utils.sam;

import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.utils.test.BaseTest;
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
