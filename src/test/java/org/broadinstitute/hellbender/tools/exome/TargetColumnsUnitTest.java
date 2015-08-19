package org.broadinstitute.hellbender.tools.exome;

import org.testng.Assert;
import org.testng.annotations.Test;

/**
 * Unit tests for {@link TargetColumns}.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public class TargetColumnsUnitTest {

    @Test
    public void testIsSpecialColumnOnSpecialNames() {
        for (final TargetColumns special : TargetColumns.values()) {
            Assert.assertTrue(TargetColumns.isTargetColumnName(special.toString()));
        }
    }

    @Test
    public void testIsSpecialColumnOnNotSpecialNames() {
        for (final TargetColumns special : TargetColumns.values()) {
            Assert.assertFalse(TargetColumns.isTargetColumnName(special.name() + "_something"));
        }
        Assert.assertFalse(TargetColumns.isTargetColumnName(""));
        Assert.assertFalse(TargetColumns.isTargetColumnName("NA12878"));
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testIsSpecialColumnOnNull() {
        TargetColumns.isTargetColumnName(null);
    }
}
