package org.broadinstitute.hellbender.tools.spark.sv;

import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Iterator;
import java.util.NoSuchElementException;

/**
 * Unit tests for SVUtils.
 */
public class SVUtilsUnitTest extends BaseTest {

    @Test(groups = "spark")
    void hashMapCapacityTest() {
        Assert.assertEquals(SVUtils.hashMapCapacity(150),201);
    }
}
