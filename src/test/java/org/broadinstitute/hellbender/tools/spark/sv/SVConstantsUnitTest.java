package org.broadinstitute.hellbender.tools.spark.sv;

import org.broadinstitute.hellbender.utils.ClassUtils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

public final class SVConstantsUnitTest extends BaseTest{
    @Test
    public void testAccess() throws Exception {
        Assert.assertFalse(ClassUtils.canMakeInstances(SVConstants.class));
    }
}
