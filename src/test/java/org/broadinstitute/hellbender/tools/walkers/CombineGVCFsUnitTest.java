package org.broadinstitute.hellbender.tools.walkers;

import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;
import java.util.stream.Collectors;

public class CombineGVCFsUnitTest {

    @DataProvider(name="breakIntermediateStopSites")
    public Object[][] getIntermediateStopSitesData() {
        return new Object[][] {
                // Note that the expected results here do not represent a final set of stop sites for the given
                // interval. Rather, they are intended to match the output of the CombineGVCFs.getIntermediateStopSites
                // method, which returns an initial set of intermediate stop sites that in some cases includes sites
                // outside the actual interval being closed, but which  are subsequently filtered out by additional
                // downstream code in CombineGVCFs.
                { new SimpleInterval("contig", 1, 1), 1, Collections.EMPTY_LIST },
                { new SimpleInterval("contig", 1, 2), 1, Arrays.asList(1) },
                { new SimpleInterval("contig", 1, 2), 2, Arrays.asList(1) },
                { new SimpleInterval("contig", 1, 10), 2, Arrays.asList(1, 3, 5, 7, 9) },
                { new SimpleInterval("contig", 1, 10), 5, Arrays.asList(4, 9) },
                { new SimpleInterval("contig", 1, 100), 25, Arrays.asList(24, 49, 74, 99) },

                { new SimpleInterval("contig", 10, 10), 2, Arrays.asList(9) },
                { new SimpleInterval("contig", 10, 10), 5, Arrays.asList(9) },

                { new SimpleInterval("contig", 10, 10), 10, Arrays.asList(9) },
                { new SimpleInterval("contig", 10, 10), 100, Collections.EMPTY_LIST },
                { new SimpleInterval("contig", 10, 20), 5, Arrays.asList(9, 14, 19) },
                { new SimpleInterval("contig", 10, 20), 10, Arrays.asList(9, 19) },
                { new SimpleInterval("contig", 10, 20), 100, Collections.EMPTY_LIST },

                { new SimpleInterval("contig", 10, 100), 25, Arrays.asList(24, 49, 74, 99) },
                { new SimpleInterval("contig", 10, 100), 50, Arrays.asList(49, 99) },
                { new SimpleInterval("contig", 10, 100), 100, Arrays.asList(99) },
                { new SimpleInterval("contig", 10, 100), 1000, Collections.EMPTY_LIST },
                { new SimpleInterval("contig", 110, 120), 100, Arrays.asList(99) }
        };
    }

    @Test(dataProvider = "breakIntermediateStopSites")
    public void testGetIntermediateStopSites(
            final SimpleInterval intervalToClose,
            final int breakBandMultiple,
            final List<Integer> expectedCloseSites)
    {
        final List<Integer> actualStopSites = new ArrayList<>(CombineGVCFs.getIntermediateStopSites(intervalToClose, breakBandMultiple));
        actualStopSites.sort(Comparator.naturalOrder());
        // validate that the resulting stop sites all result in valid single-position stop intervals
        actualStopSites.stream().forEach(stopSite -> Assert.assertNotNull(new SimpleInterval(intervalToClose.getContig(), stopSite, stopSite)));
        Assert.assertEquals(actualStopSites, expectedCloseSites);
    }

}
