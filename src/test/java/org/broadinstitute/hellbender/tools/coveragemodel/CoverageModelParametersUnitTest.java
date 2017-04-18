package org.broadinstitute.hellbender.tools.coveragemodel;

import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.annotations.Test;

/**
 * Unit tests for {@link CoverageModelParameters}
 *
 * TODO github/gatk-protected issue #843
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public class CoverageModelParametersUnitTest extends BaseTest {

    /**
     * Create a random model, write, read, and assert
     */
    @Test(enabled = false)
    public void testReadWriteModel() { }

    /**
     * Create a random model, create a random read count collection with re-arranged
     * targets, adapt, and assert
     */
    @Test(enabled = false)
    public void testAdaptModelToReadCountCollection() { }
}
