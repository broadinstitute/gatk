package org.broadinstitute.hellbender.tools.coveragemodel;

import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.List;

/**
 * Unit tests for {@link CoverageModelSparkUtils}
 *
 * TODO github/gatk-protected issue #843
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public class CoverageModelSparkUtilsUnitTest extends BaseTest {

    @Test(dataProvider = "createLinearSpaceBlocksTestData")
    public void testCreateLinearSpaceBlocks(final int length, final int numBlocks, final int minBlockSize) {
        final List<LinearlySpacedIndexBlock> blocks = CoverageModelSparkUtils.createLinearlySpacedIndexBlocks(length,
                numBlocks, minBlockSize);
        /* assert full coverage */
        final int actualNumBlocks = blocks.size();
        final LinearlySpacedIndexBlock firstBlock = blocks.get(0);
        final LinearlySpacedIndexBlock lastBlock = blocks.get(actualNumBlocks - 1);
        Assert.assertEquals(firstBlock.getBegIndex(), 0);
        Assert.assertEquals(lastBlock.getEndIndex(), length);
        for (int i = 0; i < actualNumBlocks - 1; i++) {
            Assert.assertEquals(blocks.get(i).getEndIndex(), blocks.get(i + 1).getBegIndex());
        }
        /* assert min block size */
        if (actualNumBlocks > 1) {
            for (final LinearlySpacedIndexBlock block : blocks) {
                Assert.assertTrue(block.getNumElements() >= minBlockSize);
            }
        }
        /* assert hash code */
        for (int i = 0; i < actualNumBlocks; i++) {
            Assert.assertTrue(blocks.get(i).hashCode() == i);
        }
    }

    @DataProvider(name = "createLinearSpaceBlocksTestData")
    public Object[][] getCreateLinearSpaceBlocksTestData() {
        return new Object[][] {
                {10, 3, 5}, /* length, numBLocks, minBlockSize */
                {123, 20, 4},
                {1000, 50, 1},
                {1234, 100, 20},
                {1235432, 234, 20},
                {123, 200, 1},
                {123, 200, 200}
        };
    }
}
