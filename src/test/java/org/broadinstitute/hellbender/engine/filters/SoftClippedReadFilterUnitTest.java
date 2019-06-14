package org.broadinstitute.hellbender.engine.filters;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

/**
 * Unit tests for the {@link SoftClippedReadFilter} class.
 * Created by jonn on 6/13/19.
 */
public class SoftClippedReadFilterUnitTest {

    private static final int CHR_COUNT = 1;
    private static final int CHR_START = 1;
    private static final int CHR_SIZE = 1000;
    private static final int GROUP_COUNT = 5;

    private final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeaderWithGroups(CHR_COUNT, CHR_START, CHR_SIZE, GROUP_COUNT);

    private GATKRead buildSAMRead(final String cigarString) {
        final Cigar cigar = TextCigarCodec.decode(cigarString);
        return ArtificialReadUtils.createArtificialRead(header, cigar);
    }

    @Test(dataProvider= "SoftClipRatioDataProvider")
    public void testOverclippedSoftClipRatioFilter(final String cigarString,
                                                   final double clipRatio,
                                                   final boolean expectedResult) {

        final SoftClippedReadFilter filter = new SoftClippedReadFilter();
        filter.minimumSoftClippedRatio = clipRatio;

        final GATKRead read = buildSAMRead(cigarString);
        Assert.assertEquals(filter.test(read), expectedResult, cigarString);

        filter.doInvertFilter = true;
        Assert.assertEquals(filter.test(read), !expectedResult, "Inverted case: " + cigarString);
    }

    @Test(dataProvider= "SoftClippedLeadingTrailingRatioDataProvider")
    public void testSoftClippedLeadingTrailingRatioFilter(final String cigarString,
                                                          final double clipRatio,
                                                          final boolean expectedResult) {

        final SoftClippedReadFilter filter = new SoftClippedReadFilter();
        filter.minimumLeadingTrailingSoftClippedRatio = clipRatio;

        final GATKRead read = buildSAMRead(cigarString);
        Assert.assertEquals(filter.test(read), expectedResult, cigarString);

        filter.doInvertFilter = true;
        Assert.assertEquals(filter.test(read), !expectedResult, "Inverted case: " + cigarString);
    }

    @DataProvider(name = "SoftClipRatioDataProvider")
    public Iterator<Object[]> softClipRatioDataProvider() {
        final List<Object[]> testData = new LinkedList<>();

        // ---------------------------------------
        // Null / trivial cases:
        testData.add(new Object[] { "", 0.1, false });
        testData.add(new Object[] { "10H", 0.1, false });

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

    @DataProvider(name = "SoftClippedLeadingTrailingRatioDataProvider")
    public Iterator<Object[]> softClippedLeadingTrailingRatioDataProvider() {
        final List<Object[]> testData = new LinkedList<>();

        // ---------------------------------------
        // Null / trivial cases:
        testData.add(new Object[] { "", 0.1, false });
        testData.add(new Object[] { "10H", 0.1, false });

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

}
