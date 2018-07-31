package org.broadinstitute.hellbender.utils.downsampling;

import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Collection;
import java.util.LinkedList;
import java.util.List;

public final class LevelingDownsamplerUnitTest extends GATKBaseTest {

    private static final class LevelingDownsamplerUniformStacksTest extends TestDataProvider {
        public enum DataStructure { LINKED_LIST, ARRAY_LIST }

        int targetSize;
        int numStacks;
        int stackSize;
        DataStructure dataStructure;
        int expectedSize;

        public LevelingDownsamplerUniformStacksTest(final int targetSize, final int numStacks, final int stackSize, final DataStructure dataStructure ) {
            super(LevelingDownsamplerUniformStacksTest.class);

            this.targetSize = targetSize;
            this.numStacks = numStacks;
            this.stackSize = stackSize;
            this.dataStructure = dataStructure;
            expectedSize = calculateExpectedDownsampledStackSize();

            setName(String.format("%s: targetSize=%d numStacks=%d stackSize=%d dataStructure=%s expectedSize=%d",
                    getClass().getSimpleName(), targetSize, numStacks, stackSize, dataStructure, expectedSize));
        }

        public Collection<List<Object>> createStacks() {
            final Collection<List<Object>> stacks = new ArrayList<List<Object>>();

            for ( int i = 1; i <= numStacks; i++ ) {
                final List<Object> stack = dataStructure == DataStructure.LINKED_LIST ? new LinkedList<>() : new ArrayList<>();

                for ( int j = 1; j <= stackSize; j++ ) {
                    stack.add(new Object());
                }

                stacks.add(stack);
            }

            return stacks;
        }

        private int calculateExpectedDownsampledStackSize() {
            final int numItemsToRemove = numStacks * stackSize - targetSize;

            if ( numStacks == 0 ) {
                return 0;
            }
            else if ( numItemsToRemove <= 0 ) {
                return stackSize;
            }

            return Math.max(1, stackSize - (numItemsToRemove / numStacks));
        }
    }

    @DataProvider(name = "UniformStacksDataProvider")
    public Object[][] createUniformStacksTestData() {
        for ( int targetSize = 1; targetSize <= 10000; targetSize *= 10 ) {
            for ( int numStacks = 0; numStacks <= 10; numStacks++ ) {
                for ( int stackSize = 1; stackSize <= 1000; stackSize *= 10 ) {
                    for ( final LevelingDownsamplerUniformStacksTest.DataStructure dataStructure : LevelingDownsamplerUniformStacksTest.DataStructure.values() ) {
                        new LevelingDownsamplerUniformStacksTest(targetSize, numStacks, stackSize, dataStructure);
                    }
                }
            }
        }

        return LevelingDownsamplerUniformStacksTest.getTests(LevelingDownsamplerUniformStacksTest.class);
    }

    @Test( dataProvider = "UniformStacksDataProvider" )
    public void testLevelingDownsamplerWithUniformStacks(final LevelingDownsamplerUniformStacksTest test ) {
        logger.warn("Running test: " + test);

        Utils.resetRandomGenerator();

        final Downsampler<List<Object>> downsampler = new LevelingDownsampler<>(test.targetSize);

        final Collection<List<Object>> stacks = test.createStacks();
        if (stacks.size() == 1){
            downsampler.submit(stacks.iterator().next()); //This is only done to exercise the 1 element code path in 'submit'
        } else {
            downsampler.submit(stacks);
        }


        if ( test.numStacks > 0 ) {
            Assert.assertFalse(downsampler.hasFinalizedItems());
            Assert.assertTrue(downsampler.peekFinalized() == null);
            Assert.assertTrue(downsampler.hasPendingItems());
            Assert.assertTrue(downsampler.peekPending() != null);
        }
        else {
            Assert.assertFalse(downsampler.hasFinalizedItems() || downsampler.hasPendingItems());
            Assert.assertTrue(downsampler.peekFinalized() == null && downsampler.peekPending() == null);
        }

        downsampler.signalEndOfInput();

        if ( test.numStacks > 0 ) {
            Assert.assertTrue(downsampler.hasFinalizedItems());
            Assert.assertTrue(downsampler.peekFinalized() != null);
            Assert.assertFalse(downsampler.hasPendingItems());
            Assert.assertTrue(downsampler.peekPending() == null);
        }
        else {
            Assert.assertFalse(downsampler.hasFinalizedItems() || downsampler.hasPendingItems());
            Assert.assertTrue(downsampler.peekFinalized() == null && downsampler.peekPending() == null);
        }

        final int sizeFromDownsampler = downsampler.size();
        final List<List<Object>> downsampledStacks = downsampler.consumeFinalizedItems();
        Assert.assertFalse(downsampler.hasFinalizedItems() || downsampler.hasPendingItems());
        Assert.assertTrue(downsampler.peekFinalized() == null && downsampler.peekPending() == null);

        Assert.assertEquals(downsampledStacks.size(), test.numStacks);

        int totalRemainingItems = 0;
        for ( final List<Object> stack : downsampledStacks ) {
            Assert.assertTrue(Math.abs(stack.size() - test.expectedSize) <= 1);
            totalRemainingItems += stack.size();
        }

        Assert.assertEquals(sizeFromDownsampler, totalRemainingItems);
        final int numItemsReportedDiscarded = downsampler.getNumberOfDiscardedItems();
        final int numItemsActuallyDiscarded = test.numStacks * test.stackSize - totalRemainingItems;

        Assert.assertEquals(numItemsReportedDiscarded, numItemsActuallyDiscarded);

        downsampler.resetStats();
        Assert.assertEquals(downsampler.getNumberOfDiscardedItems(), 0);

        Assert.assertTrue(totalRemainingItems <= Math.max(test.targetSize, test.numStacks));
    }
}
