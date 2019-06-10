package org.broadinstitute.hellbender.engine.filters;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;


public class OverclippedReadFilterUnitTest extends GATKBaseTest {

    private static final int CHR_COUNT = 1;
    private static final int CHR_START = 1;
    private static final int CHR_SIZE = 1000;
    private static final int GROUP_COUNT = 5;

    private final SAMFileHeader header= ArtificialReadUtils.createArtificialSamHeaderWithGroups(CHR_COUNT, CHR_START, CHR_SIZE, GROUP_COUNT);

    private GATKRead buildSAMRead(final String cigarString) {
        final Cigar cigar = TextCigarCodec.decode(cigarString);
        return ArtificialReadUtils.createArtificialRead(header, cigar);
    }

    @Test(dataProvider= "OverclippedDataProvider")
    public void testOverclippedFilter(final String cigarString,
                                      final boolean doNotRequireSoftClipsOnBothEnds,
                                      final boolean expectedResult) {

        final OverclippedReadFilter filter = new OverclippedReadFilter(30, doNotRequireSoftClipsOnBothEnds);
        final GATKRead read = buildSAMRead(cigarString);
        Assert.assertEquals(filter.test(read), expectedResult, cigarString);
    }

    @Test(dataProvider= "OverclippedSoftClipRatioDataProvider")
    public void testOverclippedSoftClipRatioFilter(final String cigarString,
                                                   final double clipRatio,
                                                   final boolean expectedResult) {

        final OverclippedReadFilter filter = new OverclippedReadFilter(10);
        filter.minimumSoftClippedRatio = clipRatio;

        final GATKRead read = buildSAMRead(cigarString);
        Assert.assertEquals(filter.test(read), expectedResult, cigarString);
    }

    @Test(dataProvider= "OverclippedLeadingTrailingSoftClipRatioDataProvider")
    public void testOverclippedLeadingTrailingSoftClipRatioFilter(final String cigarString,
                                                   final double clipRatio,
                                                   final boolean expectedResult) {

        final OverclippedReadFilter filter = new OverclippedReadFilter(10);
        filter.minimumLeadingTrailingSoftClippedRatio = clipRatio;

        final GATKRead read = buildSAMRead(cigarString);
        Assert.assertEquals(filter.test(read), expectedResult, cigarString);
    }

    @DataProvider(name = "OverclippedSoftClipRatioDataProvider")
    public Iterator<Object[]> overclippedSoftClipRatioDataProvider() {
        final List<Object[]> testData = new LinkedList<>();

        // ---------------------------------------
        // Null / trivial cases:
        testData.add(new Object[] { "", 0.1, false });
        testData.add(new Object[] { "10H", 0.1, false });

        // ---------------------------------------
        // Min bases test:
        testData.add(new Object[] { "20S5M", 0.2, false }); // 20/25 = .8
        testData.add(new Object[] { "20S20M", 0.2, true }); // 20/40 = .5

        // ---------------------------------------
        // Soft clip ratio test:

        testData.add(new Object[] { "1S1M1S17M", 0.2, false }); // 2/20 = .100
        testData.add(new Object[] { "1S1M2S17M", 0.2, false }); // 3/21 = .143
        testData.add(new Object[] { "1S1M3S17M", 0.2, false }); // 4/22 = .182
        testData.add(new Object[] { "1S1M4S17M", 0.2, true });  // 5/23 = .217
        testData.add(new Object[] { "1S1M5S17M", 0.2, true });  // 6/24 = .250
        testData.add(new Object[] { "1S1M6S17M", 0.2, true });  // 7/25 = .280

        // ---------------------------------------
        // Soft clip placement:

        testData.add(new Object[] { "101S100M", 0.5, true });
        testData.add(new Object[] { "100M101S", 0.5, true });
        testData.add(new Object[] { "25H20S10M20S10M20S10M20S10M20S10M20S25H", 0.5, true });

        return testData.iterator();
    }

    @DataProvider(name = "OverclippedLeadingTrailingSoftClipRatioDataProvider")
    public Iterator<Object[]> overclippedLeadingTrailingSoftClipRatioDataProvider() {
        final List<Object[]> testData = new LinkedList<>();

        // ---------------------------------------
        // Null / trivial cases:
        testData.add(new Object[] { "", 0.1, false });
        testData.add(new Object[] { "10H", 0.1, false });

        // ---------------------------------------
        // Min bases test:
        testData.add(new Object[] { "20S5M", 0.2, false }); // 20/25 = .8
        testData.add(new Object[] { "20S20M", 0.2, true }); // 20/40 = .5

        // ---------------------------------------
        // Soft clip ratio test:

        // Non-leading/-trailing
        testData.add(new Object[] { "1S1M1S17M", 0.2, false }); // 2/20 = .100
        testData.add(new Object[] { "1S1M2S17M", 0.2, false }); // 3/21 = .143
        testData.add(new Object[] { "1S1M3S17M", 0.2, false }); // 4/22 = .182
        testData.add(new Object[] { "1S1M4S17M", 0.2, false });  // 5/23 = .217
        testData.add(new Object[] { "1S1M5S17M", 0.2, false });  // 6/24 = .250
        testData.add(new Object[] { "1S1M6S17M", 0.2, false });  // 7/25 = .280

        // Leading:
        testData.add(new Object[] { "2S1S1S16M", 0.2, false });  // 2/20 = .100
        testData.add(new Object[] { "3S1S1S16M", 0.2, false });  // 3/21 = .143
        testData.add(new Object[] { "4S1S1S16M", 0.2, false });  // 4/22 = .182
        testData.add(new Object[] { "5S1S1S16M", 0.2, true });   // 5/23 = .217
        testData.add(new Object[] { "6S1S1S16M", 0.2, true });   // 6/24 = .250
        testData.add(new Object[] { "7S1S1S16M", 0.2, true });   // 7/25 = .280

        // Trailing:
        testData.add(new Object[] { "1M1S16M2S", 0.2, false });  // 2/20 = .100
        testData.add(new Object[] { "1M1S16M3S", 0.2, false });  // 3/21 = .143
        testData.add(new Object[] { "1M1S16M4S", 0.2, false });  // 4/22 = .182
        testData.add(new Object[] { "1M1S16M5S", 0.2, true });   // 5/23 = .217
        testData.add(new Object[] { "1M1S16M6S", 0.2, true });   // 6/24 = .250
        testData.add(new Object[] { "1M1S16M7S", 0.2, true });   // 7/25 = .280

        // ---------------------------------------
        // Soft clip placement:

        testData.add(new Object[] { "101S100M", 0.5, true });
        testData.add(new Object[] { "100M101S", 0.5, true });
        testData.add(new Object[] { "25H20S10M20S10M20S10M20S10M20S10M20S25H", 0.5, false });

        return testData.iterator();
    }

    @DataProvider(name = "OverclippedDataProvider")
    public Iterator<Object[]> overclippedDataProvider() {
        final List<Object[]> result = new LinkedList<>();

        result.add(new Object[] { "1S10M1S", false, false });
        result.add(new Object[] { "1S10X1S", false, false });
        result.add(new Object[] { "1H1S10M1S1H", false, false });
        result.add(new Object[] { "1S40M1S", false, true});
        result.add(new Object[] { "1S40X1S", false, true });
        result.add(new Object[] { "1H10M1S", false, true});
        result.add(new Object[] { "1S10M1H", false, true});

        result.add(new Object[] { "10M1S", false, true});
        result.add(new Object[] { "1S10M", false, true});

        result.add(new Object[] { "10M1S", true, false});
        result.add(new Object[] { "1S10M", true, false});

        result.add(new Object[] { "1S10M10D10M1S", false, false });
        result.add(new Object[] { "1S1M40I1S", false, true });

        result.add(new Object[] { "1S10I1S", false, false });
        result.add(new Object[] { "1S40I1S", false, true });
        result.add(new Object[] { "1S40I1S", true, true });

        result.add(new Object[] { "25S40I25M", true, true });

        //Read is too short once soft-clipping removed
        result.add(new Object[] { "25S25M", true, false });
        result.add(new Object[] { "25S25X", true, false });
        result.add(new Object[] { "25S25H", true, false });
        result.add(new Object[] { "25S25H", false, true });

        result.add(new Object[] { "25S25M25S", false, false });
        result.add(new Object[] { "25M25S", true, false });
        result.add(new Object[] { "25S25M", true, false });

        result.add(new Object[] { "25S35S", true, false });

        //Read long enough even with soft clipping removed
        result.add(new Object[] { "25S35M25S", true, true });
        result.add(new Object[] { "35M25S", true, true });
        result.add(new Object[] { "25S35M", true, true });

        return result.iterator();
    }
}
