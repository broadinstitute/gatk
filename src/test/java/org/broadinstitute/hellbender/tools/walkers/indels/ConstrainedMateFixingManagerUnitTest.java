/*
* Copyright 2012-2016 Broad Institute, Inc.
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

package org.broadinstitute.gatk.tools.walkers.indels;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.ProgressLoggerInterface;
import org.broadinstitute.gatk.utils.BaseTest;
import org.broadinstitute.gatk.utils.GenomeLocParser;
import org.broadinstitute.gatk.utils.sam.ArtificialSAMUtils;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.List;


public class ConstrainedMateFixingManagerUnitTest extends BaseTest {

    private static SAMFileHeader header;
    private static GenomeLocParser genomeLocParser;

    @BeforeClass
    public void beforeClass() {
        header = ArtificialSAMUtils.createArtificialSamHeader(3, 1, 10000);
        genomeLocParser = new GenomeLocParser(header.getSequenceDictionary());
    }

    @Test
    public void testSecondaryAlignmentsDoNotInterfere() {
        final List<GATKSAMRecord> properReads = ArtificialSAMUtils.createPair(header, "foo", 1, 10, 30, true, false);
        final GATKSAMRecord read1 = properReads.get(0);
        read1.setAlignmentStart(8); // move the read
        read1.setFlags(99);   // first in proper pair, mate negative strand

        final GATKSAMRecord read2Primary = properReads.get(1);
        read2Primary.setFlags(147);   // second in pair, mate unmapped, not primary alignment

        Assert.assertEquals(read1.getInferredInsertSize(), 21);

        final GATKSAMRecord read2NonPrimary = new GATKSAMRecord(read2Primary);
        read2NonPrimary.setFlags(393);   // second in proper pair, on reverse strand

        final ConstrainedMateFixingManager manager = new ConstrainedMateFixingManager(null, genomeLocParser, 1000, 1000, 1000);
        manager.addRead(read1, true, false);
        manager.addRead(read2NonPrimary, false, false);
        manager.addRead(read2Primary, false, false);

        Assert.assertEquals(manager.getNReadsInQueue(), 3);

        for ( final SAMRecord read : manager.getReadsInQueueForTesting() ) {
            if ( read.getFirstOfPairFlag() ) {
                Assert.assertEquals(read.getFlags(), 99);
                Assert.assertEquals(read.getInferredInsertSize(), 23);
            } else if ( read.getNotPrimaryAlignmentFlag() ) {
                Assert.assertEquals(read.getFlags(), 393);
                Assert.assertEquals(read.getInferredInsertSize(), -21);
            } else {
                Assert.assertEquals(read.getFlags(), 147);
                Assert.assertEquals(read.getInferredInsertSize(), -23);
            }
        }
    }

    @Test
    public void testSecondaryAlignmentsDoNotCauseAccidentalRemovalOfMate() {
        final List<GATKSAMRecord> properReads = ArtificialSAMUtils.createPair(header, "foo", 1, 530, 1594, true, false);
        final GATKSAMRecord read1 = properReads.get(0);
        read1.setFlags(99);   // first in proper pair, mate negative strand

        final GATKSAMRecord read2Primary = properReads.get(1);
        read2Primary.setFlags(147);   // second in pair, mate unmapped, not primary alignment
        read2Primary.setAlignmentStart(1596); // move the read

        final GATKSAMRecord read2NonPrimary = new GATKSAMRecord(read2Primary);
        read2NonPrimary.setReadName("foo");
        read2NonPrimary.setFlags(393);   // second in proper pair, on reverse strand
        read2NonPrimary.setAlignmentStart(451);
        read2NonPrimary.setMateAlignmentStart(451);

        final ConstrainedMateFixingManager manager = new ConstrainedMateFixingManager(null, genomeLocParser, 10000, 200, 10000);
        manager.addRead(read2NonPrimary, false, false);
        manager.addRead(read1, false, false);

        for ( int i = 0; i < ConstrainedMateFixingManager.EMIT_FREQUENCY; i++ )
            manager.addRead(ArtificialSAMUtils.createArtificialRead(header, "foo" + i, 0, 1500, 10), false, false);

        Assert.assertTrue(manager.forMateMatching.containsKey("foo"));
    }

    @Test
    public void testSupplementaryAlignmentsDoNotCauseBadMateFixing() {
        final List<GATKSAMRecord> properReads = ArtificialSAMUtils.createPair(header, "foo", 1, 1000, 2000, true, false);
        final GATKSAMRecord read1 = properReads.get(0);
        read1.setFlags(99);   // first in pair, negative strand

        final GATKSAMRecord read2 = properReads.get(1);
        read2.setFlags(161);   // second in pair, mate negative strand

        final GATKSAMRecord read2Supp = new GATKSAMRecord(read2);
        read2Supp.setReadName("foo");
        read2Supp.setFlags(2209);   // second in pair, mate negative strand, supplementary
        read2Supp.setAlignmentStart(100);
        read2Supp.setMateAlignmentStart(1000);

        final DummyWriter writer = new DummyWriter();
        final ConstrainedMateFixingManager manager = new ConstrainedMateFixingManager(writer, genomeLocParser, 10000, 200, 10000);
        manager.addRead(read2Supp, false, false);
        manager.addRead(read1, false, false);
        manager.addRead(read2, false, false);
        manager.close(); // "write" the reads to our dummy writer

        // check to make sure that none of the mate locations were changed, which is the problem brought to us by a user
        for ( final SAMRecord read : writer.reads ) {
            final int start = read.getAlignmentStart();
            switch (start) {
                case 100:
                    Assert.assertEquals(read.getMateAlignmentStart(), 1000);
                    break;
                case 1000:
                    Assert.assertEquals(read.getMateAlignmentStart(), 2000);
                    break;
                case 2000:
                    Assert.assertEquals(read.getMateAlignmentStart(), 1000);
                    break;
                default:
                    Assert.assertTrue(false, "We saw a read located at the wrong position");
            }
        }
    }

    private class DummyWriter implements SAMFileWriter {

        public List<SAMRecord> reads;

        public DummyWriter() { reads = new ArrayList<>(10); }

        public void addAlignment(final SAMRecord alignment) { reads.add(alignment);}

        public SAMFileHeader getFileHeader() { return null; }

        public void setProgressLogger(final ProgressLoggerInterface progress) {}

        public void close() {}
    }
}