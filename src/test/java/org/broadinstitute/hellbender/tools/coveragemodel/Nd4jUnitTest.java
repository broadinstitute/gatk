package org.broadinstitute.hellbender.tools.coveragemodel;

import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.nd4j.linalg.api.buffer.DataBuffer;
import org.nd4j.linalg.factory.Nd4j;
import org.testng.Assert;
import org.testng.annotations.Test;

/**
 * This test is to ensure that the DType is correctly set to Double
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public class Nd4jUnitTest extends BaseTest {
    @Test
    public void testDType() {
        Assert.assertTrue(Nd4j.dataType().equals(DataBuffer.Type.DOUBLE), "Data type for Nd4j must be set to" +
                " double; otherwise, coverage model EM algorithm will not function properly");
    }
}
