package org.broadinstitute.hellbender.utils.runtime;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.Main;
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

    @Test
    public void testGetToolkitName(){
        final String toolkitName = RuntimeUtils.getToolkitName(Main.class);
        // test that it's one of the valid outcomes of getToolkitName, tests in intellij don't necessarily buildAndWriteLine the correct manifest but the docker tests might\
        Assert.assertTrue(toolkitName.equals("The Genome Analysis Toolkit (GATK)") || toolkitName.equals("org.broadinstitute.hellbender"));
    }
}
