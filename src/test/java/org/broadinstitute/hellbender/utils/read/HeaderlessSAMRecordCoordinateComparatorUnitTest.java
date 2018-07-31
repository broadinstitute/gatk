package org.broadinstitute.hellbender.utils.read;

import htsjdk.samtools.*;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class HeaderlessSAMRecordCoordinateComparatorUnitTest extends GATKBaseTest {

    /**
     * Tests that the ordering produced by HeaderlessSAMRecordCoordinateComparator matches the ordering
     * produced by SAMRecordCoordinateComparator for a representative selection of reads.
     */
    @Test
    public void testComparatorOrderingMatchesHtsjdk() throws IOException {
        // This file has unmapped reads that are set to the position of their mates -- the ordering check
        // in the test below will fail if our ordering of these reads relative to the mapped reads
        // is not consistent with the definition of coordinate sorting as defined in
        // htsjdk.samtools.SAMRecordCoordinateComparator
        final String inputBam = toolsTestDir + "BQSR/CEUTrio.HiSeq.WGS.b37.ch20.1m-1m1k.NA12878.bam";
        final List<SAMRecord> originalReads = new ArrayList<>();
        final List<SAMRecord> headerlessReads = new ArrayList<>();
        SAMFileHeader header = null;

        try ( final SamReader reader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File(inputBam)) ) {
            header = reader.getFileHeader();

            for ( SAMRecord read : reader ) {
                // Clear the indexing bin so that it doesn't affect the equality checks below
                read.setFlags(read.getFlags());
                originalReads.add(read);

                final SAMRecord copy = read.deepCopy();
                copy.setHeaderStrict(null);
                // Clear the indexing bin so that it doesn't affect the equality checks below
                copy.setFlags(copy.getFlags());
                headerlessReads.add(copy);
            }
        }

        final List<SAMRecord> actualReads = new ArrayList<>(headerlessReads);
        Collections.shuffle(actualReads);
        Collections.sort(actualReads, new HeaderlessSAMRecordCoordinateComparator(header));

        final List<SAMRecord> expectedReads = new ArrayList<>(originalReads);
        Collections.shuffle(expectedReads);
        Collections.sort(expectedReads, new SAMRecordCoordinateComparator());

        Assert.assertEquals(actualReads.size(), expectedReads.size());
        for ( int i = 0; i < actualReads.size(); ++i ) {
            final SAMRecord actualRead = actualReads.get(i);
            final SAMRecord expectedRead = expectedReads.get(i);

            // Restore the header on the headerless read before the equality check
            actualRead.setHeaderStrict(header);

            Assert.assertEquals(actualRead, expectedRead, "Ordering produced by HeaderlessSAMRecordCoordinateComparator does not match the ordering produced by SAMRecordCoordinateComparator");
        }
    }
}
