package org.broadinstitute.hellbender.cmdline;

import org.broadinstitute.hellbender.Main;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Collections;

public class DeprecatedToolsRegistryTest {

    @Test
    public void testRunMissingDeprecatedTool() {
        final String missingTool = "IndelRealigner";

        final UserException e = Assert.expectThrows(
                UserException.class,
                () -> new Main().instanceMain( new String[] {missingTool} )
        );
        Assert.assertTrue(e.getMessage().contains(DeprecatedToolsRegistry.getToolDeprecationInfo(missingTool)));
    }

    @Test
    public void testRunMissingButNotRegisteredTool() {
        final String missingButNotRegisteredTool = "MadeUpToolNotInTheRegistry";

        Assert.assertNull(DeprecatedToolsRegistry.getToolDeprecationInfo(missingButNotRegisteredTool));

        final Main main = new Main();
        final UserException e = Assert.expectThrows(
                UserException.class,
                () -> main.instanceMain( new String[] {missingButNotRegisteredTool} )
        );
        Assert.assertTrue(e.getMessage().contains(main.getUnknownCommandMessage(Collections.emptySet(), missingButNotRegisteredTool)));
    }

}
