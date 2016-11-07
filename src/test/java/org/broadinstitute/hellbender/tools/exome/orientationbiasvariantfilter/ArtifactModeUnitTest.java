package org.broadinstitute.hellbender.tools.exome.orientationbiasvariantfilter;

import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

public class ArtifactModeUnitTest extends BaseTest {

    @Test
    public void testBasicCreation() {
        final ArtifactMode createdFromConstructor = new ArtifactMode('C', 'T');
        final ArtifactMode createdFromOf = ArtifactMode.of('C', 'T');

        Assert.assertEquals(createdFromConstructor, createdFromOf);
        Assert.assertEquals(createdFromConstructor.getLeft().charValue(), 'C');
        Assert.assertEquals(createdFromConstructor.getRight().charValue(), 'T');
        Assert.assertEquals(createdFromOf.getLeft().charValue(), 'C');
        Assert.assertEquals(createdFromOf.getRight().charValue(), 'T');
    }
}
