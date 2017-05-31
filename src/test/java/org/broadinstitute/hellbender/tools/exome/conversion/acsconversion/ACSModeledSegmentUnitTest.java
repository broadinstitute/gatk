package org.broadinstitute.hellbender.tools.exome.conversion.acsconversion;

import org.broadinstitute.hellbender.tools.exome.ModeledSegment;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

/**
 * Created by lichtens on 6/16/16.
 */
public class ACSModeledSegmentUnitTest extends BaseTest {
    @Test
    public void testSegMean() {
        final ACSModeledSegment tmp = new ACSModeledSegment(new SimpleInterval("1", 100, 110), ModeledSegment.NO_CALL, 10, ParamUtils.log2(1.5),
                10, .3, .1, 1, .1, 1, .1, 2);

        Assert.assertEquals(tmp.getTau(), 3.0, 1e-8);
    }
}
