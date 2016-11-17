package org.broadinstitute.hellbender.tools.spark.pathseq;

import org.broadinstitute.hellbender.tools.spark.transforms.markduplicates.MarkDuplicatesSparkUtils;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.Assert;
import org.testng.annotations.Test;

public class MarkedOpticalDuplicateReadFilterTest {

    @Test
    public void testTest() throws Exception {

        final MarkedOpticalDuplicateReadFilter filter = new MarkedOpticalDuplicateReadFilter();

        final GATKRead read_in = ArtificialReadUtils.createArtificialRead("*");
        Assert.assertEquals(filter.test(read_in),true);

        read_in.setAttribute(MarkDuplicatesSparkUtils.OPTICAL_DUPLICATE_TOTAL_ATTRIBUTE_NAME,1);
        Assert.assertEquals(filter.test(read_in),false);

        read_in.setAttribute(MarkDuplicatesSparkUtils.OPTICAL_DUPLICATE_TOTAL_ATTRIBUTE_NAME,0);
        Assert.assertEquals(filter.test(read_in),true);
    }

}