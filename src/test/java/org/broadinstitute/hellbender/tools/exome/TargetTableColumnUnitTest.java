package org.broadinstitute.hellbender.tools.exome;

import org.testng.Assert;
import org.testng.annotations.Test;

/**
 * Unit tests for {@link TargetTableColumn}.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public class TargetTableColumnUnitTest {

    @Test
    public void testIsSpecialColumnOnSpecialNames() {
        for (final TargetTableColumn special : TargetTableColumn.values()) {
            Assert.assertTrue(TargetTableColumn.isStandardTargetColumnName(special.toString()));
        }
    }

    @Test
    public void testIsSpecialColumnOnNotSpecialNames() {
        for (final TargetTableColumn special : TargetTableColumn.values()) {
            Assert.assertFalse(TargetTableColumn.isStandardTargetColumnName(special.name() + "_something"));
        }
        Assert.assertFalse(TargetTableColumn.isStandardTargetColumnName(""));
        Assert.assertFalse(TargetTableColumn.isStandardTargetColumnName("NA12878"));
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testIsSpecialColumnOnNull() {
        TargetTableColumn.isStandardTargetColumnName(null);
    }
}
