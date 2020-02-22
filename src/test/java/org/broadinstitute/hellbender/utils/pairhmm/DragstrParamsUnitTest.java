package org.broadinstitute.hellbender.utils.pairhmm;

import org.broadinstitute.hellbender.testutils.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

/**
 * Unit test for class {@link DragstrParams}.
 */
public class DragstrParamsUnitTest extends BaseTest {

    public static final String TEST_PARAMS_FILE = "org/broadinstitute/hellbender/utils/dragstr/dna_nexus_novaseq_plus0_0_params.txt";


    @Test
    public void testLoad() {
        final String filePath = ("src/test/resources/" +  TEST_PARAMS_FILE);
        final DragstrParams params = new DragstrParams(filePath);
        Assert.assertNotNull(params);
        Assert.assertEquals(params.maximumPeriod(), 8);
        Assert.assertEquals(params.maximumRepeats(), 20);
    }

    @Test
    public void testQueries() {
        final String filePath = ("src/test/resources/" +  TEST_PARAMS_FILE);
        final DragstrParams params = new DragstrParams(filePath);
        for (int i = 1; i <= params.maximumPeriod(); i++) {
            for (int j = 1; j <= params.maximumRepeats(); i++) {
                final double gop = params.gop(i, j);
                final double gcp = params.gcp(i, j);
                final double api = params.api(i, j);
                Assert.assertTrue(gop > 0.0  && !Double.isInfinite(gop));
                Assert.assertTrue(gcp > 0.0 && !Double.isInfinite(gcp));
                Assert.assertTrue(api > 0.0 && !Double.isInfinite(api));

            }
        }
    }

}
