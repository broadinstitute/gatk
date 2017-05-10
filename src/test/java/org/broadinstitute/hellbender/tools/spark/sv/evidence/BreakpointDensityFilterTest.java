package org.broadinstitute.hellbender.tools.spark.sv.evidence;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;
import org.broadinstitute.hellbender.tools.spark.utils.FlatMapGluer;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.stream.Collectors;

import static org.broadinstitute.hellbender.tools.spark.sv.evidence.BreakpointEvidence.ReadEvidence.TemplateEnd;

public class BreakpointDensityFilterTest extends BaseTest {
    private static final SAMFileHeader artificialSamHeader =
            ArtificialReadUtils.createArtificialSamHeaderWithGroups(2, 1, 1000000, 1);
    private static final ReadMetadata readMetadata = initMetadata();
    private static final PartitionCrossingChecker emptyCrossingChecker = new PartitionCrossingChecker();

    private static ReadMetadata initMetadata() {
        final ReadMetadata.PartitionBounds[] partitionBounds = new ReadMetadata.PartitionBounds[3];
        partitionBounds[0] = new ReadMetadata.PartitionBounds(0, 1, 0, 10000);
        partitionBounds[1] = new ReadMetadata.PartitionBounds(0, 10001, 0, 20000);
        partitionBounds[2] = new ReadMetadata.PartitionBounds(0, 20001, 0, 30000);
        return new ReadMetadata(new HashSet<>(), artificialSamHeader,
                                new ReadMetadata.LibraryFragmentStatistics(350, 40, 40),
                                partitionBounds, 100, 10, 30);
    }

    @DataProvider(name = "simpleEvidenceClusters")
    public Object[][] simpleEvidence() {
        final List<BreakpointEvidence> evidenceList = new ArrayList<>(5);

        final BreakpointEvidence splitEvidence1 = new BreakpointEvidence.SplitRead(ArtificialReadUtils.createArtificialRead(artificialSamHeader,
                "1", "1", 1000,
                ArtificialReadUtils.createRandomReadBases(151, false),
                ArtificialReadUtils.createRandomReadQuals(151),
                "100M50S"), readMetadata, false);
        evidenceList.add(splitEvidence1);

        final BreakpointEvidence splitEvidence2 = new BreakpointEvidence.SplitRead(ArtificialReadUtils.createArtificialRead(artificialSamHeader,
                "1", "1", 1010,
                ArtificialReadUtils.createRandomReadBases(151, false),
                ArtificialReadUtils.createRandomReadQuals(151),
                "90M60S"), readMetadata, false);
        evidenceList.add(splitEvidence2);

        final List<GATKRead> pair1 = ArtificialReadUtils.createPair(artificialSamHeader, "pair1", 151, 1250, 500000, true, false);
        evidenceList.add(new BreakpointEvidence.WeirdTemplateSize(pair1.get(0), readMetadata));

        final List<GATKRead> pair2 = ArtificialReadUtils.createPair(artificialSamHeader, "pair1", 151, 1255, 500000, true, false);
        evidenceList.add(new BreakpointEvidence.WeirdTemplateSize(pair2.get(0), readMetadata));

        final List<GATKRead> pair3 = ArtificialReadUtils.createPair(artificialSamHeader, "pair1", 151, 1350, 500000, true, false);
        evidenceList.add(new BreakpointEvidence.WeirdTemplateSize(pair3.get(0), readMetadata));

        return new Object[][] {{evidenceList}};
    }


    @Test(dataProvider = "simpleEvidenceClusters")
    public void testGetBreakpointClusters(final List<BreakpointEvidence> evidenceList) {

        BreakpointDensityFilter breakpointDensityFilter =
                new BreakpointDensityFilter(evidenceList.iterator(), readMetadata, 3, 3, emptyCrossingChecker);

        Assert.assertFalse(breakpointDensityFilter.hasEnoughOverlappers(evidenceList.get(0).getLocation()));
        Assert.assertFalse(breakpointDensityFilter.hasEnoughOverlappers(evidenceList.get(1).getLocation()));
        Assert.assertTrue(breakpointDensityFilter.hasEnoughOverlappers(evidenceList.get(2).getLocation()));
        Assert.assertTrue(breakpointDensityFilter.hasEnoughOverlappers(evidenceList.get(3).getLocation()));
        Assert.assertTrue(breakpointDensityFilter.hasEnoughOverlappers(evidenceList.get(4).getLocation()));

        Assert.assertTrue(breakpointDensityFilter.hasNext());
        Assert.assertEquals(breakpointDensityFilter.next(), evidenceList.get(2));
        Assert.assertTrue(breakpointDensityFilter.hasNext());
        Assert.assertEquals(breakpointDensityFilter.next(), evidenceList.get(3));
        Assert.assertTrue(breakpointDensityFilter.hasNext());
        Assert.assertEquals(breakpointDensityFilter.next(), evidenceList.get(4));
        Assert.assertFalse(breakpointDensityFilter.hasNext());
    }

    @Test(dataProvider = "simpleEvidenceClusters")
    public void testGetBreakpointClustersWithCoherentEvidence(final List<BreakpointEvidence> evidenceList) {

        BreakpointDensityFilter breakpointDensityFilter =
                new BreakpointDensityFilter(evidenceList.iterator(), readMetadata, 5, 3, emptyCrossingChecker);

        Assert.assertFalse(breakpointDensityFilter.hasEnoughOverlappers(evidenceList.get(0).getLocation()));
        Assert.assertFalse(breakpointDensityFilter.hasEnoughOverlappers(evidenceList.get(1).getLocation()));
        Assert.assertTrue(breakpointDensityFilter.hasEnoughOverlappers(evidenceList.get(2).getLocation()));
        Assert.assertTrue(breakpointDensityFilter.hasEnoughOverlappers(evidenceList.get(3).getLocation()));
        Assert.assertTrue(breakpointDensityFilter.hasEnoughOverlappers(evidenceList.get(4).getLocation()));

        Assert.assertTrue(breakpointDensityFilter.hasNext());
        Assert.assertEquals(breakpointDensityFilter.next(), evidenceList.get(2));
        Assert.assertTrue(breakpointDensityFilter.hasNext());
        Assert.assertEquals(breakpointDensityFilter.next(), evidenceList.get(3));
        Assert.assertTrue(breakpointDensityFilter.hasNext());
        Assert.assertEquals(breakpointDensityFilter.next(), evidenceList.get(4));
        Assert.assertFalse(breakpointDensityFilter.hasNext());
    }

    @Test
    public void testCrossPartitionEvidence() {
        // build evidence for partition 0
        final List<BreakpointEvidence> evidenceList0 = new ArrayList<>();
        // a valid cluster in the 1st fragment
        evidenceList0.add(makeEvidence(1));
        evidenceList0.add(makeEvidence(11));
        evidenceList0.add(makeEvidence(21));
        // a straggler in the 1st fragment
        evidenceList0.add(makeEvidence(131));
        // a straggler in the 2nd fragment
        evidenceList0.add(makeEvidence(401));
        // a valid cluster in the 2nd fragment
        evidenceList0.add(makeEvidence(551));
        evidenceList0.add(makeEvidence(561));
        evidenceList0.add(makeEvidence(571));
        // a straggler in the middle
        evidenceList0.add(makeEvidence(1001));
        // a valid cluster in the middle
        evidenceList0.add(makeEvidence(2001));
        evidenceList0.add(makeEvidence(2011));
        evidenceList0.add(makeEvidence(2021));
        // a straggler in the penultimate fragment
        evidenceList0.add(makeEvidence(9351));
        // a valid cluster in the penultimate fragment
        evidenceList0.add(makeEvidence(9501));
        evidenceList0.add(makeEvidence(9511));
        evidenceList0.add(makeEvidence(9521));
        // a valid cluster in the final fragment
        evidenceList0.add(makeEvidence(9751));
        evidenceList0.add(makeEvidence(9761));
        evidenceList0.add(makeEvidence(9771));
        // a straggler in the final fragment
        evidenceList0.add(makeEvidence(9901));

        final List<BreakpointEvidence> actualList0 = filterList(0, evidenceList0);
        final List<Integer> expectedStartList0 =
                Arrays.asList(1, 11, 21, 551, 561, 571, 2001, 2011, 2021, 9501, 9511, 9521, 9751, 9761, 9771, 9901);
        Assert.assertEquals(startList(actualList0), expectedStartList0);

        final List<BreakpointEvidence> actualList0Clustered = clusterList(0, actualList0);
        final List<Integer> expectedStartList0Clustered =
                Arrays.asList(1, 551, 2001, 9501, 9511, 9521, 9751, 9761, 9771, 9901);
        Assert.assertEquals(startList(actualList0Clustered), expectedStartList0Clustered);

        final List<BreakpointEvidence> evidenceList1 = new ArrayList<>(evidenceList0.size());
        evidenceList0.forEach(ev -> evidenceList1.add(makeEvidence(ev.getLocation().getStart()+10000)));
        evidenceList1.add(makeEvidence(19955)); // additional straggler that overlaps cluster in next partition
        final List<BreakpointEvidence> actualList1 = filterList(1, evidenceList1);
        final List<Integer> expectedStartList1 =
                Arrays.asList(10001, 10011, 10021, 10131, 10551, 10561, 10571, 12001, 12011, 12021, 19501, 19511, 19521, 19751, 19761, 19771, 19901, 19955);
        Assert.assertEquals(startList(actualList1), expectedStartList1);

        final List<BreakpointEvidence> actualList1Clustered = clusterList(1, actualList1);
        final List<Integer> expectedStartList1Clustered =
                Arrays.asList(10001, 10011, 10021, 10131, 10551, 10561, 10571, 12001, 19501, 19511, 19521, 19751, 19761, 19771, 19901, 19955);
        Assert.assertEquals(startList(actualList1Clustered), expectedStartList1Clustered);

        final List<BreakpointEvidence> evidenceList2 = new ArrayList<>(evidenceList0.size());
        evidenceList0.forEach(ev -> evidenceList2.add(makeEvidence(ev.getLocation().getStart()+20000)));
        final List<BreakpointEvidence> actualList2 = filterList(2, evidenceList2);
        final List<Integer> expectedStartList2 =
                Arrays.asList(20001, 20011, 20021, 20131, 20551, 20561, 20571, 22001, 22011, 22021, 29501, 29511, 29521, 29751, 29761, 29771);
        Assert.assertEquals(startList(actualList2), expectedStartList2);

        final List<BreakpointEvidence> actualList2Clustered = clusterList(2, actualList2);
        final List<Integer> expectedStartList2Clustered =
                Arrays.asList(20001, 20011, 20021, 20131, 20551, 20561, 20571, 22001, 29501);
        Assert.assertEquals(startList(actualList2Clustered), expectedStartList2Clustered);

        final List<BreakpointEvidence> allPartitions = new ArrayList<>();
        allPartitions.addAll(actualList0Clustered);
        allPartitions.addAll(actualList1Clustered);
        allPartitions.addAll(actualList2Clustered);
        BreakpointDensityFilter breakpointDensityFilter =
                new BreakpointDensityFilter(allPartitions.iterator(), readMetadata, 3, 3, emptyCrossingChecker);
        final List<BreakpointEvidence> actualAllPartitions = new ArrayList<>();
        while ( breakpointDensityFilter.hasNext() ) {
            actualAllPartitions.add(breakpointDensityFilter.next());
        }
        final List<Integer> expectedAllPartitions =
                Arrays.asList(1, 551, 2001, 9501, 9511, 9521, 9751, 9761, 9771,
                        10001, 10011, 10021, 10551, 10561, 10571, 12001, 19501, 19511, 19521, 19751, 19761, 19771, 19955,
                        20001, 20011, 20021, 20551, 20561, 20571, 22001, 29501);
        Assert.assertEquals(startList(actualAllPartitions), expectedAllPartitions);

        final List<BreakpointEvidence> actualAllPartitionsClustered = clusterList(-1, actualAllPartitions);
        final List<Integer> expectedAllPartitionsClustered =
                Arrays.asList(1, 551, 2001, 9501, 10551, 12001, 19501, 20551, 22001, 29501);
        Assert.assertEquals(startList(actualAllPartitionsClustered), expectedAllPartitionsClustered);
    }

    private List<Integer> startList( final List<BreakpointEvidence> evList ) {
        return evList.stream().map(ev -> ev.getLocation().getStart()).collect(Collectors.toList());
    }

    private List<BreakpointEvidence> clusterList( final int partitionIdx, final List<BreakpointEvidence> inputList ) {
        final int fragmentLength = readMetadata.getMaxMedianFragmentSize();
        final PartitionCrossingChecker partitionCrossingChecker = partitionIdx < 0 ? new PartitionCrossingChecker() :
                new PartitionCrossingChecker(partitionIdx,readMetadata,2*fragmentLength);
        BreakpointEvidenceClusterer breakpointEvidenceClusterer =
                new BreakpointEvidenceClusterer(fragmentLength, partitionCrossingChecker);
        final BreakpointEvidence sentinel = new BreakpointEvidence(new SVInterval(1,1,1),0,false);
        FlatMapGluer<BreakpointEvidence, BreakpointEvidence> gluer =
                new FlatMapGluer<>(breakpointEvidenceClusterer,inputList.iterator(),sentinel);
        List<BreakpointEvidence> outputList = new ArrayList<>();
        while ( gluer.hasNext() ) {
            outputList.add(gluer.next());
        }
        return outputList;
    }

    private List<BreakpointEvidence> filterList( final int partitionIdx, final List<BreakpointEvidence> inputList )  {
        final PartitionCrossingChecker crossingChecker =
                new PartitionCrossingChecker(partitionIdx, readMetadata, readMetadata.getMaxMedianFragmentSize());
        BreakpointDensityFilter breakpointDensityFilter =
                new BreakpointDensityFilter(inputList.iterator(), readMetadata, 3, 3, crossingChecker);
        final List<BreakpointEvidence> actualList = new ArrayList<>();
        while ( breakpointDensityFilter.hasNext() ) {
            actualList.add(breakpointDensityFilter.next());
        }
        return actualList;
    }

    private BreakpointEvidence.ReadEvidence makeEvidence( final int start ) {
        return new BreakpointEvidence.ReadEvidence(new SVInterval(0,start,start+100),1,"Test",TemplateEnd.UNPAIRED,false);
    }
}
