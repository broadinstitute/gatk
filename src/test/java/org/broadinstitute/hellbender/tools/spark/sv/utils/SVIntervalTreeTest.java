package org.broadinstitute.hellbender.tools.spark.sv.utils;

import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Iterator;

public class SVIntervalTreeTest extends BaseTest {
    private static SVInterval[] intervals = {
            new SVInterval(0, 100, 200),
            new SVInterval(0, 150, 250),
            new SVInterval(0, 200, 300),
            new SVInterval(0, 250, 350),
            new SVInterval(0, 300, 400),
            new SVInterval(0, 350, 450),
            new SVInterval(0, 400, 500),
            new SVInterval(0, 450, 550),
            new SVInterval(0, 500, 600),
            new SVInterval(0, 550, 650),
            new SVInterval(0, 600, 700),
            new SVInterval(0, 650, 750),
            new SVInterval(0, 700, 800),
            new SVInterval(0, 750, 850),
            new SVInterval(0, 800, 900),
            new SVInterval(0, 850, 950),
            new SVInterval(0, 900, 1000),
            new SVInterval(1, 0, 100)
    };

    private static SVIntervalTree<Integer> initTree() {
        final SVIntervalTree<Integer> tree = new SVIntervalTree<>();
        final int[] loadOrder = {1, 5, 3, 15, 11, 7, 13, 9, 14, 4, 17, 12, 0, 8, 16, 2, 6, 10};
        Assert.assertEquals(intervals.length, loadOrder.length);
        for ( final int idx : loadOrder ) {
            tree.put(intervals[idx], idx);
        }
        return tree;
    }

    @Test(groups = "sv")
    public void sizeAndClearTest() {
        final SVIntervalTree<Integer> testTree = initTree();
        Assert.assertEquals(testTree.size(), intervals.length);
        testTree.clear();
        Assert.assertEquals(testTree.size(), 0);
    }

    @Test(groups = "sv")
    public void putFindAndOrderTest() {
        final SVIntervalTree<Integer> testTree = initTree();
        for ( int idx = 0; idx != intervals.length; ++idx ) {
            final SVInterval interval = intervals[idx];
            final SVIntervalTree.Entry<Integer> entry = testTree.find(interval);
            Assert.assertEquals(entry.getInterval(), interval);
            Assert.assertEquals(entry.getValue().intValue(), idx);
        }
        final Integer sentinel = -47;
        testTree.setSentinel(sentinel);
        Assert.assertEquals(testTree.getSentinel(), sentinel);

        // putting an existing interval into the tree ought to return the old value,
        // replace the value, and leave the size unchanged
        final int someIndex = 3;
        final Integer someValue = -1;
        Assert.assertEquals(testTree.put(intervals[someIndex], someValue).intValue(), someIndex);
        Assert.assertEquals(testTree.find(intervals[someIndex]).getValue(), someValue);
        Assert.assertEquals(testTree.size(), intervals.length);

        // putting a new interval into the tree ought to return the sentinel value and increase the size by 1
        Assert.assertEquals(testTree.put(new SVInterval(-1, 0, 0), someValue), sentinel);
        Assert.assertEquals(testTree.size(), intervals.length+1);
    }

    @Test(groups = "sv")
    public void removeTest() {
        for ( int deletedIdx = 0; deletedIdx != intervals.length; ++deletedIdx ) {
            final SVIntervalTree<Integer> testTree = initTree();
            testTree.remove(intervals[deletedIdx]);
            int idx = 0;
            for ( final SVIntervalTree.Entry<Integer> entry : testTree ) {
                if ( idx == deletedIdx ) {
                    idx += 1;
                }
                Assert.assertEquals(entry.getInterval(), intervals[idx++]);
            }
        }

        // delete everything via an iterator
        final SVIntervalTree<Integer> tree = initTree();
        final Iterator<SVIntervalTree.Entry<Integer>> itr = tree.iterator();
        while ( itr.hasNext() ) {
            itr.next();
            itr.remove();
        }
        Assert.assertEquals(tree.size(), 0);
        Assert.assertFalse(tree.iterator().hasNext());
    }

    @Test(groups = "sv")
    public void findByIndexTest() {
        final SVIntervalTree<Integer> testTree = initTree();
        for ( int idx = 0; idx != intervals.length; ++idx ) {
            SVIntervalTree.Entry<Integer> entry = testTree.findByIndex(idx);
            Assert.assertEquals(entry.getInterval(), intervals[idx]);
            Assert.assertEquals(entry.getValue().intValue(), idx);
        }
    }

    @Test(groups = "sv")
    public void getIndexTest() {
        final SVIntervalTree<Integer> testTree = initTree();
        for ( int idx = 0; idx != intervals.length; ++idx ) {
            Assert.assertEquals(testTree.getIndex(intervals[idx]), idx);
            SVIntervalTree.Entry<Integer> entry = testTree.findByIndex(idx);
            Assert.assertEquals(entry.getInterval(), intervals[idx]);
            Assert.assertEquals(entry.getValue().intValue(), idx);
        }
    }

    @Test(groups = "sv")
    public void minTest() {
        final SVIntervalTree<Integer> testTree = initTree();

        Assert.assertEquals(testTree.min().getInterval(), intervals[0]);

        for ( int idx = 0; idx != intervals.length; ++idx ) {
            Assert.assertEquals(testTree.min(intervals[idx]).getInterval(), intervals[idx]);
        }
        // try an interval that's less than anything in the set
        Assert.assertEquals(testTree.min(new SVInterval(0, 0, 1)).getInterval(), intervals[0]);

        // try an interval that's greater than anything in the set
        Assert.assertNull(testTree.min(new SVInterval(2, 0, 1)));

        // try an empty interval
        Assert.assertEquals(testTree.min(new SVInterval(0, 201, 201)).getInterval(), intervals[3]);
    }

    @Test(groups = "sv")
    public void minOverlapperTest() {
        final SVIntervalTree<Integer> testTree = initTree();
        Assert.assertEquals(testTree.minOverlapper(intervals[0]).getInterval(), intervals[0]);
        final int lastIdx = intervals.length - 1;
        for ( int idx = 1; idx != lastIdx; ++idx ) {
            Assert.assertEquals(testTree.minOverlapper(intervals[idx]).getInterval(), intervals[idx-1]);
        }
        Assert.assertEquals(testTree.minOverlapper(intervals[lastIdx]).getInterval(), intervals[lastIdx]);
        // try an interval that's less than anything in the set
        Assert.assertNull(testTree.minOverlapper(new SVInterval(0, 0, 1)));

        // try an interval that's greater than anything in the set
        Assert.assertNull(testTree.minOverlapper(new SVInterval(2, 0, 1)));

        // try an empty interval
        Assert.assertEquals(testTree.minOverlapper(new SVInterval(0, 201, 201)).getInterval(), intervals[1]);
    }

    @Test(groups = "sv")
    public void maxTest() {
        final SVIntervalTree<Integer> testTree = initTree();

        final int lastIdx = intervals.length - 1;
        Assert.assertEquals(testTree.max().getInterval(), intervals[lastIdx]);

        for ( int idx = 0; idx != intervals.length; ++idx ) {
            Assert.assertEquals(testTree.max(intervals[idx]).getInterval(), intervals[idx]);
        }
        // try an interval that's less than anything in the set
        Assert.assertNull(testTree.max(new SVInterval(0, 0, 1)));

        // try an interval that's greater than anything in the set
        Assert.assertEquals(testTree.max(new SVInterval(2, 0, 1)).getInterval(), intervals[lastIdx]);

        // try an empty interval
        Assert.assertEquals(testTree.max(new SVInterval(0, 201, 201)).getInterval(), intervals[2]);
    }

    @Test(groups = "sv")
    public void iteratorTest() {
        final SVIntervalTree<Integer> testTree = initTree();
        final Iterator<SVIntervalTree.Entry<Integer>> itr1 = testTree.iterator();
        int idx = 0;
        while ( itr1.hasNext() ) {
            Assert.assertEquals(itr1.next().getInterval(), intervals[idx++]);
        }
        Assert.assertEquals(idx, intervals.length);

        idx = 3;
        final Iterator<SVIntervalTree.Entry<Integer>> itr2 = testTree.iterator(intervals[idx]);
        while ( itr2.hasNext() ) {
            Assert.assertEquals(itr2.next().getInterval(), intervals[idx++]);
        }

        // try an empty tree
        Assert.assertFalse(new SVIntervalTree<Integer>().iterator().hasNext());
    }

    @Test(groups = "sv")
    public void overlappersTest() {
        final SVIntervalTree<Integer> testTree = initTree();
        final Iterator<SVIntervalTree.Entry<Integer>> itr1 = testTree.overlappers(new SVInterval(0, 240, 260));
        int idx = 1;
        while ( itr1.hasNext() ) {
            Assert.assertEquals(itr1.next().getInterval(), intervals[idx++]);
        }
        Assert.assertEquals(idx, 4);

        // test overlappers in the presence of an element in the middle that isn't an overlapper
        final SVIntervalTree<Integer> tree = new SVIntervalTree<>();
        tree.put(new SVInterval(0, 100, 1000), 0);
        tree.put(new SVInterval(0, 200, 1000), 0);
        tree.put(new SVInterval(0, 300, 1000), 0);
        tree.put(new SVInterval(0, 350, 450), 0);
        tree.put(new SVInterval(0, 400, 1000), 0);
        tree.put(new SVInterval(0, 500, 1000), 0);
        final Iterator<SVIntervalTree.Entry<Integer>> itr2 = tree.overlappers(new SVInterval(0, 500, 600));
        int start = 100;
        while ( itr2.hasNext() ) {
            Assert.assertEquals(itr2.next().getInterval().getStart(), start);
            start += 100;
        }
        Assert.assertEquals(start, 600);
    }

    @Test(groups = "sv")
    public void reverseIteratorTest() {
        final SVIntervalTree<Integer> testTree = initTree();
        final Iterator<SVIntervalTree.Entry<Integer>> itr1 = testTree.reverseIterator();
        int idx = intervals.length;
        while ( itr1.hasNext() ) {
            Assert.assertEquals(itr1.next().getInterval(), intervals[--idx]);
        }
        Assert.assertEquals(idx, 0);

        idx = 3;
        final Iterator<SVIntervalTree.Entry<Integer>> itr2 = testTree.reverseIterator(intervals[idx]);
        while ( itr2.hasNext() ) {
            Assert.assertEquals(itr2.next().getInterval(), intervals[idx--]);
        }
        Assert.assertEquals(idx, -1);

        // try an empty tree
        Assert.assertFalse(new SVIntervalTree<Integer>().reverseIterator().hasNext());
    }
}
