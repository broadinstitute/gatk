package org.broadinstitute.hellbender.cmdline;

import org.broadinstitute.hellbender.utils.ClassUtils;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

public final class StandardArgumentDefinitionsUnitTest extends GATKBaseTest {
    @Test
    public void testAccess() throws Exception {
        Assert.assertFalse(ClassUtils.canMakeInstances(StandardArgumentDefinitions.class));
    }
}
