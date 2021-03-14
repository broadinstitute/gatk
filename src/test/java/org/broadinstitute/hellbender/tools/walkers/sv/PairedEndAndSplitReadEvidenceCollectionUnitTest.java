package org.broadinstitute.hellbender.tools.walkers.sv;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.tools.sv.SplitReadEvidence;
import org.broadinstitute.hellbender.utils.codecs.SplitReadEvidenceCodec;
import org.broadinstitute.hellbender.utils.io.FeatureOutputStream;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.mockito.Mockito;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.HashSet;
import java.util.Map;
import java.util.PriorityQueue;
import java.util.function.Function;
import java.util.stream.Collectors;

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

    // Workaround for use of generic class with Mockito
    private class SplitReadFeatureOutputStream extends FeatureOutputStream<SplitReadEvidence> {
        public SplitReadFeatureOutputStream() {
            super(new GATKPath(""), new SplitReadEvidenceCodec(), SplitReadEvidenceCodec::encode,
                    new SAMSequenceDictionary(), 4);
        }
    }

    @Test
    public void testCountSplitRead() throws Exception {
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader(2, 1, 10000);
        final GATKRead rightClip = ArtificialReadUtils.createArtificialRead(header, "rightClip", 0, 1000, ArtificialReadUtils.createRandomReadBases(151, false),
                ArtificialReadUtils.createRandomReadQuals(151), "100M51S");

        final FeatureOutputStream<SplitReadEvidence> mockSrWriter = Mockito.mock(SplitReadFeatureOutputStream.class);

        PairedEndAndSplitReadEvidenceCollection tool = new PairedEndAndSplitReadEvidenceCollection();
        final PriorityQueue<PairedEndAndSplitReadEvidenceCollection.SplitPos> splitCounts = new PriorityQueue<>(new PairedEndAndSplitReadEvidenceCollection.SplitPosComparator());
        tool.sampleName = "sample";

        tool.countSplitRead(rightClip, splitCounts, mockSrWriter);
        Map<PairedEndAndSplitReadEvidenceCollection.SplitPos, Long> counts = splitCounts.stream().collect(Collectors.groupingBy(Function.identity(), Collectors.counting()));
        Assert.assertEquals(counts.get(new PairedEndAndSplitReadEvidenceCollection.SplitPos(1100, PairedEndAndSplitReadEvidenceCollection.POSITION.RIGHT)).intValue(), 1);
        Mockito.verifyZeroInteractions(mockSrWriter);

        final GATKRead rightClip2 = ArtificialReadUtils.createArtificialRead(header, "rightClip2", 0, 1050, ArtificialReadUtils.createRandomReadBases(151, false),
                ArtificialReadUtils.createRandomReadQuals(151), "50M101S");
        tool.countSplitRead(rightClip2, splitCounts, mockSrWriter);
        counts = splitCounts.stream().collect(Collectors.groupingBy(Function.identity(), Collectors.counting()));
        Assert.assertEquals(counts.get(new PairedEndAndSplitReadEvidenceCollection.SplitPos(1100, PairedEndAndSplitReadEvidenceCollection.POSITION.RIGHT)).intValue(), 2);
        Mockito.verifyZeroInteractions(mockSrWriter);

        final GATKRead leftClip = ArtificialReadUtils.createArtificialRead(header, "leftClip", 0, 1100, ArtificialReadUtils.createRandomReadBases(151, false),
                ArtificialReadUtils.createRandomReadQuals(151), "20S131M");
        tool.countSplitRead(leftClip, splitCounts, mockSrWriter);
        counts = splitCounts.stream().collect(Collectors.groupingBy(Function.identity(), Collectors.counting()));
        Assert.assertEquals(counts.get(new PairedEndAndSplitReadEvidenceCollection.SplitPos(1100, PairedEndAndSplitReadEvidenceCollection.POSITION.RIGHT)).intValue(), 2);
        Assert.assertEquals(counts.get(new PairedEndAndSplitReadEvidenceCollection.SplitPos(1100, PairedEndAndSplitReadEvidenceCollection.POSITION.LEFT)).intValue(), 1);
        Mockito.verifyZeroInteractions(mockSrWriter);

        final GATKRead leftClipDownstream = ArtificialReadUtils.createArtificialRead(header, "leftClipDownstream", 0, 1600, ArtificialReadUtils.createRandomReadBases(151, false),
                ArtificialReadUtils.createRandomReadQuals(151), "20S131M");
        tool.countSplitRead(leftClipDownstream, splitCounts, mockSrWriter);
        counts = splitCounts.stream().collect(Collectors.groupingBy(Function.identity(), Collectors.counting()));
        Assert.assertFalse(counts.containsKey(new PairedEndAndSplitReadEvidenceCollection.SplitPos(1100, PairedEndAndSplitReadEvidenceCollection.POSITION.RIGHT)));
        Assert.assertFalse(counts.containsKey(new PairedEndAndSplitReadEvidenceCollection.SplitPos(1100, PairedEndAndSplitReadEvidenceCollection.POSITION.LEFT)));
        Assert.assertEquals(counts.get(new PairedEndAndSplitReadEvidenceCollection.SplitPos(1600, PairedEndAndSplitReadEvidenceCollection.POSITION.LEFT)).intValue(), 1);
        final SplitReadEvidence splitRead1 = new SplitReadEvidence("sample", "1", 1100, 1, false);
        final SplitReadEvidence splitRead2 = new SplitReadEvidence("sample", "1", 1100, 2, true);
        Mockito.verify(mockSrWriter).add(splitRead1);
        Mockito.verify(mockSrWriter).add(splitRead2);
        Mockito.verifyNoMoreInteractions(mockSrWriter);

    }

}