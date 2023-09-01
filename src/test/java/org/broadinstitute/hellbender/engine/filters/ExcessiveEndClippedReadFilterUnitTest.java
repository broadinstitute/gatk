package org.broadinstitute.hellbender.engine.filters;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;


public class ExcessiveEndClippedReadFilterUnitTest extends GATKBaseTest {

    private static final int CHR_COUNT = 1;
    private static final int CHR_START = 1;
    private static final int CHR_SIZE = 1000;
    private static final int GROUP_COUNT = 5;

    private final SAMFileHeader header= ArtificialReadUtils.createArtificialSamHeaderWithGroups(CHR_COUNT, CHR_START, CHR_SIZE, GROUP_COUNT);

    private GATKRead buildSAMRead(final String cigarString) {
        final Cigar cigar = TextCigarCodec.decode(cigarString);
        return ArtificialReadUtils.createArtificialRead(header, cigar);
    }

    @Test(dataProvider= "provideForTestExcessiveEndClipFilter")
    public void testExcessiveEndClipFilter(final int maxClippedBases,
                                           final String cigarString,
                                           final boolean expectedResult) {

        final ExcessiveEndClippedReadFilter filter = new ExcessiveEndClippedReadFilter(maxClippedBases);
        final GATKRead read = buildSAMRead(cigarString);
        Assert.assertEquals(filter.test(read), expectedResult, cigarString);
    }

    @DataProvider(name = "provideForTestExcessiveEndClipFilter")
    public Object[][] overclippedDataProvider() {
        return new Object[][] {
                // Failing test cases:
                { 50, "1S10M", true },
                { 50, "10M1S", true },
                { 50, "49S10M", true },
                { 50, "50S10M", true },
                { 50, "10M49S", true },
                { 50, "10M50S", true },
                { 50, "49S10M1S", true },
                { 50, "50S10M1S", true },
                { 50, "1S10M49S", true },
                { 50, "1S10M50S", true },
                { 50, "49H10M1S", true },
                { 50, "50H10M1S", true },
                { 50, "1H10M49S", true },
                { 50, "1H10M50S", true },
                { 50, "49H10M49H", true },
                { 50, "50H10M50H", true },
                { 50, "40H9S10M", true },
                { 50, "40H10S10M", true },
                { 50, "10M40H9S", true },
                { 50, "10M40H10S", true },

                // Passing test cases:
                { 50, "10M51S", false },
                { 50, "51S10M", false },
                { 50, "51S10M51S", false },

                { 50, "10M51H", false },
                { 50, "51H10M", false },
                { 50, "51H10M51H", false },

                { 50, "10M26S26H", false },
                { 50, "26H26S10M", false },
                { 50, "26H26S10M26S26H", false },
        };
    }
}
