package org.broadinstitute.hellbender.tools.spark.sv.evidence;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;
import org.broadinstitute.hellbender.tools.spark.sv.utils.StrandedInterval;
import org.broadinstitute.hellbender.tools.spark.utils.IntHistogram;
import org.broadinstitute.hellbender.utils.IntHistogramTest;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

import static org.testng.Assert.assertEquals;
import static org.testng.Assert.assertTrue;

public class EvidenceTargetLinkClustererTest extends GATKBaseTest {
    private static final SAMFileHeader artificialSamHeader =
            ArtificialReadUtils.createArtificialSamHeaderWithGroups(2, 1, 1000000, 1);
    private static final ReadMetadata readMetadata = initMetadata();


    // todo: simplify; this was copied from BreakpointDensityFilterTest and we don't need all this prob
    private static ReadMetadata initMetadata() {
        final ReadMetadata.PartitionBounds[] partitionBounds = new ReadMetadata.PartitionBounds[3];
        partitionBounds[0] = new ReadMetadata.PartitionBounds(0, 1, 0, 10000, 9999);
        partitionBounds[1] = new ReadMetadata.PartitionBounds(0, 10001, 0, 20000, 9999);
        partitionBounds[2] = new ReadMetadata.PartitionBounds(0, 20001, 0, 30000, 9999);
        return new ReadMetadata(Collections.emptySet(), artificialSamHeader,
                new LibraryStatistics(new IntHistogram.CDF(IntHistogramTest.genLogNormalSample(350, 40, 10000)),
                        60000000000L, 600000000L, 1200000000000L, 3000000000L),
                partitionBounds, 100, 10, 30);
    }

    @DataProvider(name = "evidence")
    public Object[][] createTestData() {

        final List<BreakpointEvidence> evidenceList = new ArrayList<>();
        final List<GATKRead> pair1 = ArtificialReadUtils.createPair(artificialSamHeader, "pair1", 151, 1250, 500000, true, false);
        evidenceList.add(new BreakpointEvidence.WeirdTemplateSize(pair1.get(0), readMetadata));

        final List<GATKRead> pair2 = ArtificialReadUtils.createPair(artificialSamHeader, "pair2", 151, 1275, 600000, true, false);
        evidenceList.add(new BreakpointEvidence.WeirdTemplateSize(pair2.get(0), readMetadata));

        final List<GATKRead> pair3 = ArtificialReadUtils.createPair(artificialSamHeader, "pair3", 151, 1300, 500025, true, false);
        evidenceList.add(new BreakpointEvidence.WeirdTemplateSize(pair3.get(0), readMetadata));

        final List<EvidenceTargetLink> linkList = new ArrayList<>();
        linkList.add(new EvidenceTargetLink(
                new StrandedInterval(new SVInterval(readMetadata.getContigID(pair1.get(0).getContig()),
                        1275 + pair2.get(0).getLength(),
                        1275 + readMetadata.getFragmentLengthStatistics(pair2.get(0).getReadGroup()).getMaxNonOutlierFragmentSize()),
                true),
                new StrandedInterval(new SVInterval(readMetadata.getContigID(pair1.get(1).getContig()),
                        600000 - readMetadata.getFragmentLengthStatistics(pair1.get(0).getReadGroup()).getMaxNonOutlierFragmentSize() + 151,
                        600000 + BreakpointEvidence.DiscordantReadPairEvidence.MATE_ALIGNMENT_LENGTH_UNCERTAINTY),
                false),
                0, 1, Collections.singleton("foo"), new HashSet<>()));

        linkList.add(new EvidenceTargetLink(
                new StrandedInterval(new SVInterval(readMetadata.getContigID(pair1.get(0).getContig()),
                        1300 + pair3.get(0).getLength(),
                        1250 + readMetadata.getFragmentLengthStatistics(pair2.get(0).getReadGroup()).getMaxNonOutlierFragmentSize()),
                true),
                new StrandedInterval(new SVInterval(readMetadata.getContigID(pair1.get(1).getContig()),
                        500025 - readMetadata.getFragmentLengthStatistics(pair2.get(0).getReadGroup()).getMaxNonOutlierFragmentSize() + 151,
                        500000 + BreakpointEvidence.DiscordantReadPairEvidence.MATE_ALIGNMENT_LENGTH_UNCERTAINTY),
                false),
                0, 2, Collections.singleton("bar"), new HashSet<>()));

        final List<Object[]> tests = new ArrayList<>();

        tests.add(new Object[] {evidenceList.iterator(), linkList.iterator()});

        return tests.toArray(new Object[][]{});
    }


    @Test(dataProvider = "evidence", groups = "sv")
    public void testClusterEvidence( final Iterator<BreakpointEvidence> evidenceIterator,
                                     final Iterator<EvidenceTargetLink> expectedResults) throws Exception {
        final EvidenceTargetLinkClusterer clusterer = new EvidenceTargetLinkClusterer(readMetadata, 0);
        final Iterator<EvidenceTargetLink> results = clusterer.cluster(evidenceIterator);

        while (expectedResults.hasNext()) {
            final EvidenceTargetLink nextExpected =  expectedResults.next();
            assertTrue(results.hasNext());
            final EvidenceTargetLink nextActual =  results.next();
            assertEquals(nextActual, nextExpected);
        }
        assertTrue(!results.hasNext());
    }

}