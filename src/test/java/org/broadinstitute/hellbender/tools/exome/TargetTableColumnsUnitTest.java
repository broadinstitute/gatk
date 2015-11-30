package org.broadinstitute.hellbender.tools.exome;

import org.testng.Assert;
import org.testng.annotations.Test;

/**
 * Unit tests for {@link TargetTableColumns}.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public class TargetTableColumnsUnitTest {

    @Test
    public void testIsSpecialColumnOnSpecialNames() {
        for (final TargetTableColumns special : TargetTableColumns.values()) {
            Assert.assertTrue(TargetTableColumns.isStandardTargetColumnName(special.toString()));
        }
    }

    @Test
    public void testIsSpecialColumnOnNotSpecialNames() {
        for (final TargetTableColumns special : TargetTableColumns.values()) {
            Assert.assertFalse(TargetTableColumns.isStandardTargetColumnName(special.name() + "_something"));
        }
        Assert.assertFalse(TargetTableColumns.isStandardTargetColumnName(""));
        Assert.assertFalse(TargetTableColumns.isStandardTargetColumnName("NA12878"));
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testIsSpecialColumnOnNull() {
        TargetTableColumns.isStandardTargetColumnName(null);
    }
}
