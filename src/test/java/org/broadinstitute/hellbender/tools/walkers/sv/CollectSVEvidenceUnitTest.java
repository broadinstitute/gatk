package org.broadinstitute.hellbender.tools.walkers.sv;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.tools.sv.SiteDepth;
import org.broadinstitute.hellbender.tools.sv.SplitReadEvidence;
import org.broadinstitute.hellbender.tools.walkers.sv.CollectSVEvidence.*;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.codecs.SplitReadEvidenceCodec;
import org.broadinstitute.hellbender.utils.io.FeatureOutputStream;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.mockito.Mockito;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;

public class CollectSVEvidenceUnitTest extends GATKBaseTest {
    @Test
    public void testLocusComparator() {
        final SAMFileHeader hdr = ArtificialReadUtils.createArtificialSamHeader();
        final SAMSequenceDictionary dict = hdr.getSequenceDictionary();
        Assert.assertTrue(dict.getSequences().size() >= 3, "need 3 sequences");
        final LocusComparator lComp = new LocusComparatorImpl(dict);
        final String curContig = dict.getSequence(1).getContig();
        final SimpleInterval interval =
                new SimpleInterval(curContig, 1001, 2000);
        final String prevContig = dict.getSequence(0).getContig();
        Assert.assertEquals(lComp.compareLocus(prevContig, 1500, interval), -1);
        final String nextContig = dict.getSequence(2).getContig();
        Assert.assertEquals(lComp.compareLocus(nextContig, 1500, interval), 1);
        Assert.assertEquals(lComp.compareLocus(curContig, 1000, interval), -1);
        Assert.assertEquals(lComp.compareLocus(curContig, 1001, interval), 0);
        Assert.assertEquals(lComp.compareLocus(curContig, 2000, interval), 0);
        Assert.assertEquals(lComp.compareLocus(curContig, 2001, interval), 1);
    }

    @Test
    public void testEffectOfCigarAndMinQOnAlleleCounting() {
        final SAMFileHeader hdr = ArtificialReadUtils.createArtificialSamHeader();
        final LocusComparator lComp = new LocusComparatorImpl(hdr.getSequenceDictionary());
        final GATKRead read = ArtificialReadUtils.createArtificialRead(hdr, "30M5D30M5I30M");
        final int readQ = read.getBaseQuality(0);
        Assert.assertEquals(read.getBase(0), (byte)'A');
        final String contig = read.getContig();
        final int start = read.getStart();
        final String sample = "smpl";
        final List<SiteDepth> sites = new ArrayList<>(5);
        sites.add(new SiteDepth(contig, start - 1, sample, 0, 0, 0, 0));
        sites.add(new SiteDepth(contig, start, sample, 0, 0, 0, 0));
        sites.add(new SiteDepth(contig, start+30, sample, 0, 0, 0, 0));
        sites.add(new SiteDepth(contig, start+35, sample, 0, 0, 0, 0));
        sites.add(new SiteDepth(contig, start+95, sample, 0, 0, 0, 0));
        AlleleCounter.walkReadMatches(read, readQ+1, sites, lComp);
        for ( final SiteDepth sd : sites ) {
            // no counts because quality too low
            Assert.assertEquals(sd.getDepth(0), 0);
        }
        AlleleCounter.walkReadMatches(read, readQ, sites, lComp);
        Assert.assertEquals(sites.get(0).getDepth(0), 0); // upstream, so no count
        Assert.assertEquals(sites.get(1).getDepth(0), 1); // in first 30M, 1 count
        Assert.assertEquals(sites.get(2).getDepth(0), 0); // in 5D, so no count
        Assert.assertEquals(sites.get(3).getDepth(0), 1); // in second 30M, 1 count
        Assert.assertEquals(sites.get(4).getDepth(0), 0); // downstream, so no count
    }

    @Test
    public void testGetReportableDiscordantReadPair() {
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

        final CollectSVEvidence tool = new CollectSVEvidence();

        DiscordantRead reportableDiscordantReadPair =
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

    private void assertDiscordantReadInfo(final DiscordantRead reportableDiscordantReadPair, final GATKRead gatkRead) {
        Assert.assertEquals(reportableDiscordantReadPair.getName(), gatkRead.getName());
        Assert.assertEquals(reportableDiscordantReadPair.getContig(), gatkRead.getContig());
        Assert.assertEquals(reportableDiscordantReadPair.getStart(), gatkRead.getStart());
        Assert.assertEquals(reportableDiscordantReadPair.getMateContig(), gatkRead.getMateContig());
        Assert.assertEquals(reportableDiscordantReadPair.getMateStart(), gatkRead.getMateStart());
    }

    // Workaround for use of generic class with Mockito
    private static class SplitReadFeatureOutputStream extends FeatureOutputStream<SplitReadEvidence> {
        public SplitReadFeatureOutputStream() {
            super(new GATKPath(""), new SplitReadEvidenceCodec().getTabixFormat(),
                    SplitReadEvidenceCodec::encode, new SAMSequenceDictionary(), 4);
        }
    }

    @Test
    public void testCountSplitRead() {
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader(2, 1, 10000);
        final GATKRead rightClip = ArtificialReadUtils.createArtificialRead(header, "rightClip", 0, 1000, ArtificialReadUtils.createRandomReadBases(151, false),
                ArtificialReadUtils.createRandomReadQuals(151), "100M51S");

        final SplitReadFeatureOutputStream mockSrWriter = Mockito.mock(SplitReadFeatureOutputStream.class);

        CollectSVEvidence tool = new CollectSVEvidence();
        final PriorityQueue<SplitPos> splitCounts = new PriorityQueue<>(new SplitPosComparator());
        tool.sampleName = "sample";

        tool.countSplitRead(rightClip, splitCounts, mockSrWriter);
        Map<SplitPos, Long> counts = splitCounts.stream().collect(Collectors.groupingBy(Function.identity(), Collectors.counting()));
        Assert.assertEquals(counts.get(new SplitPos(1100, POSITION.RIGHT)).intValue(), 1);
        Mockito.verifyZeroInteractions(mockSrWriter);

        final GATKRead rightClip2 = ArtificialReadUtils.createArtificialRead(header, "rightClip2", 0, 1050, ArtificialReadUtils.createRandomReadBases(151, false),
                ArtificialReadUtils.createRandomReadQuals(151), "50M101S");
        tool.countSplitRead(rightClip2, splitCounts, mockSrWriter);
        counts = splitCounts.stream().collect(Collectors.groupingBy(Function.identity(), Collectors.counting()));
        Assert.assertEquals(counts.get(new SplitPos(1100, POSITION.RIGHT)).intValue(), 2);
        Mockito.verifyZeroInteractions(mockSrWriter);

        final GATKRead leftClip = ArtificialReadUtils.createArtificialRead(header, "leftClip", 0, 1100, ArtificialReadUtils.createRandomReadBases(151, false),
                ArtificialReadUtils.createRandomReadQuals(151), "20S131M");
        tool.countSplitRead(leftClip, splitCounts, mockSrWriter);
        counts = splitCounts.stream().collect(Collectors.groupingBy(Function.identity(), Collectors.counting()));
        Assert.assertEquals(counts.get(new SplitPos(1100, POSITION.RIGHT)).intValue(), 2);
        Assert.assertEquals(counts.get(new SplitPos(1100, POSITION.LEFT)).intValue(), 1);
        Mockito.verifyZeroInteractions(mockSrWriter);

        final GATKRead leftClipDownstream = ArtificialReadUtils.createArtificialRead(header, "leftClipDownstream", 0, 1600, ArtificialReadUtils.createRandomReadBases(151, false),
                ArtificialReadUtils.createRandomReadQuals(151), "20S131M");
        tool.countSplitRead(leftClipDownstream, splitCounts, mockSrWriter);
        counts = splitCounts.stream().collect(Collectors.groupingBy(Function.identity(), Collectors.counting()));
        Assert.assertFalse(counts.containsKey(new SplitPos(1100, POSITION.RIGHT)));
        Assert.assertFalse(counts.containsKey(new SplitPos(1100, POSITION.LEFT)));
        Assert.assertEquals(counts.get(new SplitPos(1600, POSITION.LEFT)).intValue(), 1);
        Mockito.verify(mockSrWriter).write(new SplitReadEvidence(tool.sampleName, "1", 1100, 1, false));
        Mockito.verify(mockSrWriter).write(new SplitReadEvidence(tool.sampleName, "1", 1100, 2, true));
        Mockito.verifyNoMoreInteractions(mockSrWriter);
    }
}