package org.broadinstitute.hellbender.utils.runtime;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

public final class RuntimeUtilsUnitTest extends GATKBaseTest {
    @Test
    public void testWhichExists() {
        Assert.assertNotNull(RuntimeUtils.which("ls"), "Unable to locate ls");
    }

    @Test
    public void testWhichNotExists() {
        Assert.assertNull(RuntimeUtils.which("does_not_exist"), "Found nonexistent binary: does_not_exist");
    }
}
