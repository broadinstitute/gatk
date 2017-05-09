package org.broadinstitute.hellbender.engine.filters;

import org.broadinstitute.hellbender.tools.spark.transforms.markduplicates.MarkDuplicatesSparkUtils;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.Assert;
import org.testng.annotations.Test;

public class NotOpticalDuplicateReadFilterTest {

    @Test
    public void testOpticalDuplicates() throws Exception {

        final NotOpticalDuplicateReadFilter filter = new NotOpticalDuplicateReadFilter();
        final GATKRead read_in = ArtificialReadUtils.createArtificialRead("*");

        //No tag set
        Assert.assertEquals(filter.test(read_in), true);

        //Is optical duplicate
        read_in.setAttribute(MarkDuplicatesSparkUtils.OPTICAL_DUPLICATE_TOTAL_ATTRIBUTE_NAME, 1);
        Assert.assertEquals(filter.test(read_in), false);

        //Is optical duplicate
        read_in.setAttribute(MarkDuplicatesSparkUtils.OPTICAL_DUPLICATE_TOTAL_ATTRIBUTE_NAME, 2);
        Assert.assertEquals(filter.test(read_in), false);

        //Not optical duplicate
        read_in.setAttribute(MarkDuplicatesSparkUtils.OPTICAL_DUPLICATE_TOTAL_ATTRIBUTE_NAME, 0);
        Assert.assertEquals(filter.test(read_in), true);
    }

}