package org.broadinstitute.hellbender.utils.read;

import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.annotations.Test;

import java.util.Iterator;

import static org.testng.Assert.*;

public final class ArtificialReadUtilsUnitTest extends BaseTest {


    @Test
    public void basicReadIteratorTest() {
        Iterator<GATKRead> iter = ArtificialReadUtils.mappedReadIterator(1, 100, 100);
        int count = 0;
        while (iter.hasNext()) {
            GATKRead rec = iter.next();
            count++;
        }
        assertEquals(count, 100 * 100);
    }

    @Test
    public void tenPerChromosome() {
        ArtificialReadQueryIterator iter = ArtificialReadUtils.mappedReadIterator(1, 100, 10);
        int count = 0;
        while (iter.hasNext()) {
            GATKRead rec = iter.next();

            assertEquals(Math.round(count / 10), ReadUtils.getReferenceIndex(rec, iter.getHeader()));
            count++;
        }
        assertEquals(count, 100 * 10);
    }

    @Test
    public void onePerChromosome() {
        ArtificialReadQueryIterator iter = ArtificialReadUtils.mappedReadIterator(1, 100, 1);
        int count = 0;
        while (iter.hasNext()) {
            GATKRead rec = iter.next();

            assertEquals(count, ReadUtils.getReferenceIndex(rec, iter.getHeader()));
            count++;
        }
        assertEquals(count, 100);
    }

    @Test
    public void basicUnmappedIteratorTest() {
        Iterator<GATKRead> iter = ArtificialReadUtils.mappedAndUnmappedReadIterator(1, 100, 100, 1000);
        int count = 0;
        for (int x = 0; x < (100* 100); x++ ) {
            if (!iter.hasNext()) {
                fail ("we didn't get the expected number of reads");
            }
            GATKRead rec = iter.next();
            assertFalse(rec.isUnmapped());
            count++;
        }
        assertEquals(100 * 100, count);

        // now we should have 1000 unmapped reads
        count = 0;
        while (iter.hasNext()) {
            GATKRead rec = iter.next();
            assertTrue(rec.isUnmapped());
            count++;
        }
        assertEquals(count, 1000);
    }


}
