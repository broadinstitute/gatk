package org.broadinstitute.hellbender.tools.spark.sv;

import org.testng.Assert;
import org.testng.annotations.Test;

/**
 * Unit test for {@link InsertSizeDistributionShape}.
 */
public class InsertSizeDistributionShapeUnitTest {

    @Test
    public void testUniqueAliases() {
        for (final InsertSizeDistributionShape shape1 : InsertSizeDistributionShape.values()) {
            for (final InsertSizeDistributionShape shape2 : InsertSizeDistributionShape.values()) {
                if (shape1 != shape2) {
                    for (final String alias1 : shape1.aliases()) {
                        for (final String alias2 : shape2.aliases()) {
                            Assert.assertNotEquals(alias1.toLowerCase(), alias2.toLowerCase());
                        }
                    }
                }
            }
        }
    }

    @Test
    public void testAliasesAndNameEncodeEachShape() {
        for (final InsertSizeDistributionShape shape : InsertSizeDistributionShape.values()) {
            Assert.assertSame(InsertSizeDistributionShape.decode(shape.name()), shape);
            for (final String alias : shape.aliases()) {
                Assert.assertSame(InsertSizeDistributionShape.decode(alias), shape);
            }
        }
    }

    @Test
    public void testDecodeOnGarbageReturnsNull() {
        Assert.assertNull(InsertSizeDistributionShape.decode("Garbage"));
    }
}
