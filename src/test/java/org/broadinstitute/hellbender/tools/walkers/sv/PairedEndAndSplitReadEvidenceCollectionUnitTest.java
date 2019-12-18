package org.broadinstitute.hellbender.tools.walkers.sv;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.HashSet;

public class PairedEndAndSplitReadEvidenceCollectionUnitTest extends GATKBaseTest {

    @Test
    public void testGetReportableDiscordantReadPair() throws Exception {
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader(2, 1, 10000);

        // read pair on different contigs
        final GATKRead discRPDiffContigsFirst = ArtificialReadUtils.createArtificialRead(header, "discRP1", 0, 1000, 150);
        discRPDiffContigsFirst.setMatePosition(header.getSequence(1).getSequenceName(), 9000);
        final GATKRead discRPDiffContigsSecond = ArtificialReadUtils.createArtificialRead(header, "discRP1", 1, 9000, 150);
        discRPDiffContigsSecond.setMatePosition(header.getSequence(0).getSequenceName(), 1000);

        // read pair on same contig
        final GATKRead discRPSameContigFirst = ArtificialReadUtils.createArtificialRead(header, "discRP2", 1, 500, 150);
        discRPSameContigFirst.setMatePosition(header.getSequence(1).getSequenceName(), 5000);
        final GATKRead discRPSameContigSecond = ArtificialReadUtils.createArtificialRead(header, "discRP2", 1, 5000, 150);
        discRPSameContigSecond.setMatePosition(header.getSequence(1).getSequenceName(), 500);

        // read pair on same contig, same start position
        final GATKRead discRPSameContigSameStartFirst = ArtificialReadUtils.createArtificialRead(header, "discRP3", 1, 4000, 150);
        discRPSameContigSameStartFirst.setMatePosition(header.getSequence(1).getSequenceName(), 4000);
        final GATKRead discRPSameContigSameStartSecond = ArtificialReadUtils.createArtificialRead(header, "discRP3", 1, 4000, 150);
        discRPSameContigSameStartSecond.setMatePosition(header.getSequence(1).getSequenceName(), 4000);

        final HashSet<String> observedDiscordantNamesAtThisLocus = new HashSet<>();

        final PairedEndAndSplitReadEvidenceCollection tool = new PairedEndAndSplitReadEvidenceCollection();

        PairedEndAndSplitReadEvidenceCollection.DiscordantRead reportableDiscordantReadPair =
                tool.getReportableDiscordantReadPair(discRPDiffContigsFirst, observedDiscordantNamesAtThisLocus, header.getSequenceDictionary());
        Assert.assertNotNull(reportableDiscordantReadPair);
        assertDiscordantReadInfo(reportableDiscordantReadPair, discRPDiffContigsFirst);

        reportableDiscordantReadPair =
                tool.getReportableDiscordantReadPair(discRPDiffContigsSecond, observedDiscordantNamesAtThisLocus, header.getSequenceDictionary());
        Assert.assertNull(reportableDiscordantReadPair);

        reportableDiscordantReadPair =
                tool.getReportableDiscordantReadPair(discRPSameContigFirst, observedDiscordantNamesAtThisLocus, header.getSequenceDictionary());
        assertDiscordantReadInfo(reportableDiscordantReadPair, discRPSameContigFirst);

        reportableDiscordantReadPair =
                tool.getReportableDiscordantReadPair(discRPSameContigSecond, observedDiscordantNamesAtThisLocus, header.getSequenceDictionary());
        Assert.assertNull(reportableDiscordantReadPair);

        reportableDiscordantReadPair =
                tool.getReportableDiscordantReadPair(discRPSameContigSameStartFirst, observedDiscordantNamesAtThisLocus, header.getSequenceDictionary());
        assertDiscordantReadInfo(reportableDiscordantReadPair, discRPSameContigSameStartFirst);
        Assert.assertTrue(observedDiscordantNamesAtThisLocus.contains(discRPSameContigSameStartFirst.getName()));

        reportableDiscordantReadPair =
                tool.getReportableDiscordantReadPair(discRPSameContigSameStartSecond, observedDiscordantNamesAtThisLocus, header.getSequenceDictionary());
        Assert.assertNull(reportableDiscordantReadPair);
        Assert.assertTrue(observedDiscordantNamesAtThisLocus.isEmpty());

    }

    private void assertDiscordantReadInfo(final PairedEndAndSplitReadEvidenceCollection.DiscordantRead reportableDiscordantReadPair, final GATKRead gatkRead) {
        Assert.assertEquals(reportableDiscordantReadPair.getName(), gatkRead.getName());
        Assert.assertEquals(reportableDiscordantReadPair.getContig(), gatkRead.getContig());
        Assert.assertEquals(reportableDiscordantReadPair.getStart(), gatkRead.getStart());
        Assert.assertEquals(reportableDiscordantReadPair.getMateContig(), gatkRead.getMateContig());
        Assert.assertEquals(reportableDiscordantReadPair.getMateStart(), gatkRead.getMateStart());
    }

    @Test
    public void testCountSplitRead() throws Exception {

    }
}