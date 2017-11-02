package org.broadinstitute.hellbender.tools.copynumber.utils;

import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.utils.IndexRange;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.List;

/**
 * Unit tests for {@link CachedBinarySearchIntervalList}
 */
public class CachedBinarySearchIntervalListTest {

    @DataProvider(name = "correctTestData")
    public Object[][] testData() {
        final List<SimpleInterval> sameChromosomeIntervalList = new ArrayList<>();
        sameChromosomeIntervalList.add(new SimpleInterval("10", 500, 1000));
        sameChromosomeIntervalList.add(new SimpleInterval("10", 1500, 2000));
        sameChromosomeIntervalList.add(new SimpleInterval("10", 2500, 3000));
        sameChromosomeIntervalList.add(new SimpleInterval("10", 5000, 6000));

        final SimpleInterval singleIntervalIntersectionCase = new SimpleInterval("10", 1400, 1600);
        final IndexRange singleIntervalIntersectionExpectedRange = new IndexRange(1, 2);

        final SimpleInterval multipleIntervalIntersectionCase = new SimpleInterval("10", 1900, 2600);
        final IndexRange multipleIntervalIntersectionExpectedRange = new IndexRange(1, 3);

        final SimpleInterval firstIntervalIntersectionCase = new SimpleInterval("10", 400, 600);
        final IndexRange firstIntervalIntersectionExpectedRange = new IndexRange(0, 1);

        final SimpleInterval queryLocationInsideIntervalCase = new SimpleInterval("10", 1700, 1800);
        final IndexRange queryLocationInsideIntervalExpectedRange = new IndexRange(1, 2);

        final SimpleInterval queryLocationSurroundsIntervalCase = new SimpleInterval("10", 400, 1100);
        final IndexRange queryLocationSurroundsIntervalExpectedRange = new IndexRange(0, 1);

        final SimpleInterval lastIntervalIntersectionCase = new SimpleInterval("10", 5900, 6100);
        final IndexRange lastIntervalIntersectionExpectedRange = new IndexRange(3, 4);

        final SimpleInterval singlePointLocationCase = new SimpleInterval("10", 1750, 1750);
        final IndexRange singlePointLocationExpectedRange = new IndexRange(1, 2);

        final SimpleInterval noIntersectionCase = new SimpleInterval("10", 2100, 2200);
        final IndexRange noIntersectionExpectedRange = new IndexRange(2, 2);

        final List<SimpleInterval> multipleChromosomeIntervalList = new ArrayList<>();
        multipleChromosomeIntervalList.add(new SimpleInterval("10", 500, 1000));
        multipleChromosomeIntervalList.add(new SimpleInterval("20", 1500, 2000));

        final SimpleInterval multipleChromosomesCase = new SimpleInterval("20", 1400, 1600);
        final IndexRange multipleChromosomesExpectedRange = new IndexRange(1, 2);

        return new Object[][] {
                {sameChromosomeIntervalList, singleIntervalIntersectionCase, singleIntervalIntersectionExpectedRange},
                {sameChromosomeIntervalList, multipleIntervalIntersectionCase, multipleIntervalIntersectionExpectedRange},
                {sameChromosomeIntervalList, firstIntervalIntersectionCase, firstIntervalIntersectionExpectedRange},
                {sameChromosomeIntervalList, queryLocationInsideIntervalCase, queryLocationInsideIntervalExpectedRange},
                {sameChromosomeIntervalList, queryLocationSurroundsIntervalCase, queryLocationSurroundsIntervalExpectedRange},
                {sameChromosomeIntervalList, lastIntervalIntersectionCase, lastIntervalIntersectionExpectedRange},
                {sameChromosomeIntervalList, singlePointLocationCase, singlePointLocationExpectedRange},
                {sameChromosomeIntervalList, noIntersectionCase, noIntersectionExpectedRange},
                {multipleChromosomeIntervalList, multipleChromosomesCase, multipleChromosomesExpectedRange}
        };
    }

    @Test(dataProvider = "correctTestData")
    public void test(final List<Locatable> intervalList, final Locatable location, final IndexRange expectedIndexRange) {
        final CachedBinarySearchIntervalList<Locatable> cachedBinarySearchIntervalList = new CachedBinarySearchIntervalList<>(intervalList);
        final IndexRange actualIndexRange = cachedBinarySearchIntervalList.findIntersectionRange(location);

        Assert.assertEquals(actualIndexRange, expectedIndexRange);
    }
}