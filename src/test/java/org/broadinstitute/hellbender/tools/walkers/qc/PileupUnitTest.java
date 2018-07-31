package org.broadinstitute.hellbender.tools.walkers.qc;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.pileup.PileupElement;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;

/**
 * Unit tests for {@link Pileup}.
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public final class PileupUnitTest extends GATKBaseTest {

    @Test
    public void testInsertLengthOutput() throws Exception {
        // create two artifitial reads
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader();
        final SimpleInterval loc = new SimpleInterval("1:1");
        final GATKRead read1 = ArtificialReadUtils.createArtificialRead(header, "read1", 0, 1, 10);
        read1.setFragmentLength(100);
        final GATKRead read2 = ArtificialReadUtils.createArtificialRead(header, "read2", 0, 1, 10);
        read2.setFragmentLength(50);
        // generate the pileups
        final ReadPileup pileup = new ReadPileup(loc, Arrays.asList(read1, read2), 1);
        // test one read
        Assert.assertEquals(Pileup.insertLengthOutput(pileup.makeFilteredPileup(
                p -> p.getRead().getName().equals("read1"))), "100");
        Assert.assertEquals(Pileup.insertLengthOutput(pileup.makeFilteredPileup(
                p -> p.getRead().getName().equals("read2"))), "50");
        // test two reads
        Assert.assertEquals(Pileup.insertLengthOutput(pileup), "100,50");
        // test an empty pileup
        Assert.assertEquals(Pileup.insertLengthOutput(new ReadPileup(loc)), "");
    }

    @Test
    public void testCreateVerboseOutput() throws Exception {
        // the header
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader();
        final SimpleInterval loc = new SimpleInterval("1:2");
        // create one read/pileup element with deletion with the following verbose string
        final String read1String = "read1@1@10@20";
        final GATKRead read1 = ArtificialReadUtils.createArtificialRead(header, "read1", 0, 1, 10);
        read1.setMappingQuality(20);
        read1.setCigar("1M1D9M");
        final PileupElement pe1 = new PileupElement(read1, 1, read1.getCigar().getCigarElement(1), 1, 0);
        // create a second one without it with the following verbose string
        final String read2String = "read2@1@50@10";
        final GATKRead read2 = ArtificialReadUtils.createArtificialRead(header, "read2", 0, 2, 50);
        read2.setMappingQuality(10);
        read2.setCigar("50M");
        final PileupElement pe2 = PileupElement.createPileupForReadAndOffset(read2, 1);

        // generate the pileups
        final ReadPileup pileup = new ReadPileup(loc, Arrays.asList(pe1, pe2));

        // test one read
        Assert.assertEquals(Pileup.createVerboseOutput(pileup.makeFilteredPileup(
                p -> p.getRead().getName().equals("read1"))), "1 "+read1String);
        Assert.assertEquals(Pileup.createVerboseOutput(pileup.makeFilteredPileup(
                p -> p.getRead().getName().equals("read2"))), "0 "+read2String);
        // test two reads
        Assert.assertEquals(Pileup.createVerboseOutput(pileup), "1 "+read1String+","+read2String);
        // test an empty pileup
        Assert.assertEquals(Pileup.createVerboseOutput(new ReadPileup(loc)), "0 ");
    }
}
