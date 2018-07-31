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


    @Test(dataProvider= "OverclippedDataProvider")
    public void testOverclippedFilter(final String cigarString, boolean doNotRequireSoftClipsOnBothEnds, final boolean expectedResult) {

        final OverclippedReadFilter filter = new OverclippedReadFilter(30, false);
        filter.doNotRequireSoftClipsOnBothEnds = doNotRequireSoftClipsOnBothEnds;
        final GATKRead read = buildSAMRead(cigarString);
        Assert.assertEquals(filter.test(read), expectedResult, cigarString);
    }

    private GATKRead buildSAMRead(final String cigarString) {
        final Cigar cigar = TextCigarCodec.decode(cigarString);
        return ArtificialReadUtils.createArtificialRead(header, cigar);
    }

    @DataProvider(name= "OverclippedDataProvider")
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
