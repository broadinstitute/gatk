package org.broadinstitute.hellbender.utils.iterators;

import htsjdk.samtools.*;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class SamReaderQueryingIteratorUnitTest extends GATKBaseTest {

    @DataProvider(name = "SamReaderQueryingIteratorTestData")
    public Object[][] samReaderQueryingIteratorTestData() {
        // This bam has only mapped reads
        final File mappedBam = new File(publicTestDir + "org/broadinstitute/hellbender/engine/reads_data_source_test1.bam");

        // This bam has mapped reads from various contigs, plus a few unmapped reads with no mapped mate
        final File unmappedBam = new File(publicTestDir + "org/broadinstitute/hellbender/engine/reads_data_source_test1_with_unmapped.bam");

        // This is a snippet of the CEUTrio.HiSeq.WGS.b37.NA12878 bam from large, with mapped reads
        // from chromosome 20 (with one mapped read having an unmapped mate), plus several unmapped
        // reads with no mapped mate.
        final File ceuSnippet = new File(publicTestDir + "org/broadinstitute/hellbender/engine/CEUTrio.HiSeq.WGS.b37.NA12878.snippet_with_unmapped.bam");

        return new Object[][] {
                // One interval, no unmapped
                { unmappedBam, Arrays.asList(new SimpleInterval("1", 200, 1000)), false, Arrays.asList("a", "b", "c", "d") },
                // One interval, with unmapped
                { unmappedBam, Arrays.asList(new SimpleInterval("1", 200, 1000)), true, Arrays.asList("a", "b", "c", "d", "u1", "u2", "u3", "u4", "u5") },
                // Multiple intervals, no unmapped
                { unmappedBam, Arrays.asList(new SimpleInterval("1", 200, 1000), new SimpleInterval("2", 500, 700), new SimpleInterval("4", 700, 701)), false, Arrays.asList("a", "b", "c", "d", "f", "g", "h", "k") },
                // Multiple intervals, with unmapped
                { unmappedBam, Arrays.asList(new SimpleInterval("1", 200, 1000), new SimpleInterval("2", 500, 700), new SimpleInterval("4", 700, 701)), true, Arrays.asList("a", "b", "c", "d", "f", "g", "h", "k", "u1", "u2", "u3", "u4", "u5") },
                // Interval with no overlapping reads, no unmapped
                { unmappedBam, Arrays.asList(new SimpleInterval("1", 3000, 4000)), false, Collections.<String>emptyList() },
                // Interval with no overlapping reads, with unmapped
                { unmappedBam, Arrays.asList(new SimpleInterval("1", 3000, 4000)), true, Arrays.asList("u1", "u2", "u3", "u4", "u5") },
                // Interval with no overlapping reads, with unmapped, but no unmapped reads in bam
                { mappedBam, Arrays.asList(new SimpleInterval("1", 3000, 4000)), true, Collections.<String>emptyList() },
                // Interval with overlapping reads, with unmapped, but no unmapped reads in bam
                { mappedBam, Arrays.asList(new SimpleInterval("1", 200, 1000)), true, Arrays.asList("a", "b", "c", "d") },
                // Null intervals, with unmapped
                { unmappedBam, null, true, Arrays.asList("u1", "u2", "u3", "u4", "u5") },
                // Empty intervals, with unmapped
                { unmappedBam, Collections.<SimpleInterval>emptyList(), true, Arrays.asList("u1", "u2", "u3", "u4", "u5") },
                // Null intervals, with unmapped, but no unmapped reads in bam
                { mappedBam, null, true, Collections.<String>emptyList() },
                // Empty intervals, with unmapped, but no unmapped reads in bam
                { mappedBam, Collections.<SimpleInterval>emptyList(), true, Collections.<String>emptyList() },
                // Null intervals, no unmapped
                { unmappedBam, null, false, Collections.<String>emptyList() },
                // Empty intervals, no unmapped
                { unmappedBam, Collections.<SimpleInterval>emptyList(), false, Collections.<String>emptyList() },
                // Interval containing mapped read with unmapped mate, no unmapped
                { ceuSnippet, Arrays.asList(new SimpleInterval("20", 10000011, 10000013)), false, Arrays.asList("a", "b", "c", "d", "e", "f", "f")},
                // Interval containing mapped read with unmapped mate, with unmapped
                { ceuSnippet, Arrays.asList(new SimpleInterval("20", 10000011, 10000013)), true, Arrays.asList("a", "b", "c", "d", "e", "f", "f", "g", "h", "h", "i", "i")},
                // Interval not containing mapped read with unmapped mate, no unmapped
                { ceuSnippet, Arrays.asList(new SimpleInterval("20", 10000009, 10000011)), false, Arrays.asList("a", "b", "c", "d", "e")},
                // Interval not containing mapped read with unmapped mate, with unmapped
                { ceuSnippet, Arrays.asList(new SimpleInterval("20", 10000009, 10000011)), true, Arrays.asList("a", "b", "c", "d", "e", "g", "h", "h", "i", "i")}
        };
    }

    @Test(dataProvider = "SamReaderQueryingIteratorTestData")
    public void testSamReaderQueryingIterator( final File inputBam, final List<SimpleInterval> queryIntervals, final boolean queryUnmapped, final List<String> expectedReadNames ) throws IOException {
        try ( final SamReader reader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(inputBam) ) {
            final SamReaderQueryingIterator iter = new SamReaderQueryingIterator(reader, queryIntervals, queryUnmapped);

            int readCount = 0;
            while ( iter.hasNext() ) {
                final SAMRecord read = iter.next();
                Assert.assertTrue(readCount < expectedReadNames.size(), "Too many reads returned from iterator");
                Assert.assertEquals(read.getReadName(), expectedReadNames.get(readCount), "Wrong read name for read returned from iterator");
                ++readCount;
            }

            Assert.assertEquals(readCount, expectedReadNames.size(), "Wrong number of reads returned from iterator");
        }
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testThrowOnNullReader() {
        final SamReaderQueryingIterator iter = new SamReaderQueryingIterator(null, Arrays.asList(new SimpleInterval("1", 1, 1)), false);
    }
}
