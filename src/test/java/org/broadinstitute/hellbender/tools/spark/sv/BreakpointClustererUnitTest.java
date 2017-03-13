package org.broadinstitute.hellbender.tools.spark.sv;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

public class BreakpointClustererUnitTest extends BaseTest {

    @DataProvider(name = "simpleEvidenceClusters")
    public Object[][] simpleEvidence() {
        final SAMFileHeader artificialSamHeader = ArtificialReadUtils.createArtificialSamHeaderWithGroups(2, 1, 1000000, 1);

        final List<BreakpointEvidence> evidenceList = new ArrayList<>(5);

        final ReadMetadata readMetadata = new ReadMetadata(new HashSet<Integer>(), artificialSamHeader,
                new ReadMetadata.ReadGroupFragmentStatistics(350, 40, 40),
                1, 100, 10, 30);

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

        return new Object[][] {{readMetadata, evidenceList}};
    }


    @Test(dataProvider = "simpleEvidenceClusters")
    public void testGetBreakpointClusters(final ReadMetadata readMetadata, final List<BreakpointEvidence> evidenceList) {

        BreakpointClusterer breakpointClusterer = new BreakpointClusterer(readMetadata, 3, 3, evidenceList.iterator());

        Assert.assertFalse(breakpointClusterer.hasEnoughOverlappers(evidenceList.get(0).getLocation()));
        Assert.assertFalse(breakpointClusterer.hasEnoughOverlappers(evidenceList.get(1).getLocation()));
        Assert.assertTrue(breakpointClusterer.hasEnoughOverlappers(evidenceList.get(2).getLocation()));
        Assert.assertTrue(breakpointClusterer.hasEnoughOverlappers(evidenceList.get(3).getLocation()));
        Assert.assertTrue(breakpointClusterer.hasEnoughOverlappers(evidenceList.get(4).getLocation()));

        Assert.assertTrue(breakpointClusterer.hasNext());
        Assert.assertEquals(breakpointClusterer.next(), evidenceList.get(2));
        Assert.assertTrue(breakpointClusterer.hasNext());
        Assert.assertEquals(breakpointClusterer.next(), evidenceList.get(3));
        Assert.assertTrue(breakpointClusterer.hasNext());
        Assert.assertEquals(breakpointClusterer.next(), evidenceList.get(4));
        Assert.assertFalse(breakpointClusterer.hasNext());
    }

    @Test(dataProvider = "simpleEvidenceClusters")
    public void testGetBreakpointClustersWithCoherentEvidence(final ReadMetadata readMetadata, final List<BreakpointEvidence> evidenceList) {

        BreakpointClusterer breakpointClusterer = new BreakpointClusterer(readMetadata, 5, 3, evidenceList.iterator());

        Assert.assertFalse(breakpointClusterer.hasEnoughOverlappers(evidenceList.get(0).getLocation()));
        Assert.assertFalse(breakpointClusterer.hasEnoughOverlappers(evidenceList.get(1).getLocation()));
        Assert.assertTrue(breakpointClusterer.hasEnoughOverlappers(evidenceList.get(2).getLocation()));
        Assert.assertTrue(breakpointClusterer.hasEnoughOverlappers(evidenceList.get(3).getLocation()));
        Assert.assertTrue(breakpointClusterer.hasEnoughOverlappers(evidenceList.get(4).getLocation()));

        Assert.assertTrue(breakpointClusterer.hasNext());
        Assert.assertEquals(breakpointClusterer.next(), evidenceList.get(2));
        Assert.assertTrue(breakpointClusterer.hasNext());
        Assert.assertEquals(breakpointClusterer.next(), evidenceList.get(3));
        Assert.assertTrue(breakpointClusterer.hasNext());
        Assert.assertEquals(breakpointClusterer.next(), evidenceList.get(4));
        Assert.assertFalse(breakpointClusterer.hasNext());
    }

}