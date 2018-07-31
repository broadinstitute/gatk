package org.broadinstitute.hellbender.tools.spark.sv.evidence;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.utils.IntHistogram;
import org.broadinstitute.hellbender.utils.IntHistogramTest;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class KSWindowFinderUnitTest extends GATKBaseTest {
    @Test(groups = "sv")
    public void testContigCrossing() {
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeaderWithGroups(3, 1, 10000000, 1);
        final SAMReadGroupRecord readGroup = header.getReadGroups().get(0);
        final String groupName = readGroup.getReadGroupId();
        final List<GATKRead> reads = new ArrayList<>(4);
        reads.addAll(ArtificialReadUtils.createPair(header,"template1",151,0,1,250,true,false));
        reads.addAll(ArtificialReadUtils.createPair(header,"template2",151,1,1,250,true,false));
        for ( final GATKRead read : reads ) {
            read.setMappingQuality(60);
            read.setReadGroup(groupName);
        }
        final StructuralVariationDiscoveryArgumentCollection.FindBreakpointEvidenceSparkArgumentCollection params =
                new StructuralVariationDiscoveryArgumentCollection.FindBreakpointEvidenceSparkArgumentCollection();
        final SVReadFilter filter = new SVReadFilter(params);
        final LibraryStatistics stats =
                new LibraryStatistics(IntHistogramTest.genLogNormalSample(400, 175, 10000).getCDF(),
                        60000000000L, 600000000L, 1200000000000L, 3000000000L);
        final Set<Integer> crossContigIgnoreSet = new HashSet<>(3);
        crossContigIgnoreSet.add(2);
        final ReadMetadata readMetadata = new ReadMetadata(crossContigIgnoreSet, header, stats, null, 4L, 4L, 1);
        final KSWindowFinder finder = new KSWindowFinder(readMetadata, filter);
        final List<BreakpointEvidence> evList = new ArrayList<>();
        for ( final GATKRead read : reads ) {
            finder.testReadAndGatherEvidence(read, evList);
        }
        final IntHistogram[] histoPair = finder.getLibraryToHistoPairMap().get(readGroup.getLibrary());
        Assert.assertEquals(histoPair[1].getTotalObservations(), 0L);
    }
}
