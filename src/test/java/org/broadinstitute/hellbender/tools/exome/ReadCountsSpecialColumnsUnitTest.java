package org.broadinstitute.hellbender.tools.exome;

import org.testng.Assert;
import org.testng.annotations.Test;

/**
 * Unit tests for {@link ReadCountsSpecialColumns}.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public class ReadCountsSpecialColumnsUnitTest {

    @Test
    public void testIsSpecialColumnOnSpecialNames() {
        for (final ReadCountsSpecialColumns special : ReadCountsSpecialColumns.values()) {
            Assert.assertTrue(ReadCountsSpecialColumns.isSpecialColumnName(special.name()));
        }
    }

    @Test
    public void testIsSpecialColumnOnNotSpecialNames() {
        for (final ReadCountsSpecialColumns special : ReadCountsSpecialColumns.values()) {
            Assert.assertFalse(ReadCountsSpecialColumns.isSpecialColumnName(special.name() + "_something"));
        }
        Assert.assertFalse(ReadCountsSpecialColumns.isSpecialColumnName(""));
        Assert.assertFalse(ReadCountsSpecialColumns.isSpecialColumnName("NA12878"));
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testIsSpecialColumnOnNull() {
        ReadCountsSpecialColumns.isSpecialColumnName(null);
    }
}
