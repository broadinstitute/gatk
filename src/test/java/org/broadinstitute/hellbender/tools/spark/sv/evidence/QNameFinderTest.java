package org.broadinstitute.hellbender.tools.spark.sv.evidence;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVIntervalTree;
import org.broadinstitute.hellbender.utils.IntHistogramTest;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;

public class QNameFinderTest extends GATKBaseTest {
    private final static LibraryStatistics LIBRARY_STATISTICS =
            new LibraryStatistics(IntHistogramTest.genLogNormalSample(400, 175, 10000).getCDF(),
                    60000000000L, 600000000L, 1200000000000L, 3000000000L);

    @Test
    public void testHighDepthRegionFiltering() throws Exception {
        final StructuralVariationDiscoveryArgumentCollection.FindBreakpointEvidenceSparkArgumentCollection params = new StructuralVariationDiscoveryArgumentCollection.FindBreakpointEvidenceSparkArgumentCollection();
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeaderWithGroups(3, 1, 10000000, 1);
        final Set<Integer> crossContigIgnoreSet = new HashSet<>();
        final ReadMetadata readMetadata = new ReadMetadata(crossContigIgnoreSet, header, LIBRARY_STATISTICS, null, 2L, 2L, 1);
        final ArrayList<SVInterval> intervals = new ArrayList<>(2);
        SVInterval interval1 = new SVInterval(0, 10378, 12002);
        SVInterval interval2 = new SVInterval(0, 115732072, 115733072);
        intervals.add(interval1);
        intervals.add(interval2);
        final GATKRead read1 = ArtificialReadUtils.createArtificialRead(header, "read1", 0, 11074,
                ArtificialReadUtils.createRandomReadBases(151, false),
                ArtificialReadUtils.createRandomReadQuals(151),
                "151M");
        read1.setIsReverseStrand(true);
        final GATKRead read2 = ArtificialReadUtils.createArtificialRead(header, "read2", 0, 11063,
                ArtificialReadUtils.createRandomReadBases(151, false),
                ArtificialReadUtils.createRandomReadQuals(151),
                "99M52S");

        final SVIntervalTree<SVInterval> highDepthIntervals = new SVIntervalTree<>();
        final SVInterval highDepthInterval1 = new SVInterval(0, 11010, 11590);
        final SVInterval highDepthInterval2 = new SVInterval(0, 115732072, 115733072);
        highDepthIntervals.put(highDepthInterval1, highDepthInterval1);
        highDepthIntervals.put(highDepthInterval2, highDepthInterval2);

        final QNameFinder qNameFinder = new QNameFinder(readMetadata, intervals, new SVReadFilter(params), highDepthIntervals);

        Iterator<QNameAndInterval> read1Result = qNameFinder.apply(read1);
        Assert.assertTrue(! read1Result.hasNext());

        Iterator<QNameAndInterval> read2Result = qNameFinder.apply(read2);
        Assert.assertTrue(! read2Result.hasNext());

    }

}