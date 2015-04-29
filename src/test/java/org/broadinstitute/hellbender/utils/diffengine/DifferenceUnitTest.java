package org.broadinstitute.hellbender.utils.diffengine;

import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

/**
 * Basic unit test for DifferableReaders in reduced reads
 */
public final class DifferenceUnitTest extends BaseTest {
    // --------------------------------------------------------------------------------
    //
    // testing routines
    //
    // --------------------------------------------------------------------------------

    private class DifferenceTest extends TestDataProvider {
        public DiffElement tree1, tree2;
        public String difference;

        private DifferenceTest(String tree1, String tree2, String difference) {
            this(DiffElement.fromString(tree1), DiffElement.fromString(tree2), difference);
        }

        private DifferenceTest(DiffElement tree1, DiffElement tree2, String difference) {
            super(DifferenceTest.class);
            this.tree1 = tree1;
            this.tree2 = tree2;
            this.difference = difference;
        }

        public String toString() {
            return String.format("tree1=%s tree2=%s diff=%s",
                    tree1 == null ? "null" : tree1.toOneLineString(),
                    tree2 == null ? "null" : tree2.toOneLineString(),
                    difference);
        }
    }

    @DataProvider(name = "data")
    public Object[][] createTrees() {
        new DifferenceTest("A=X", "A=Y", "A:1:X!=Y");
        new DifferenceTest("A=Y", "A=X", "A:1:Y!=X");
        new DifferenceTest(DiffElement.fromString("A=X"), null, "A:1:X!=MISSING");
        new DifferenceTest(null, DiffElement.fromString("A=X"), "A:1:MISSING!=X");
        return DifferenceTest.getTests(DifferenceTest.class);
    }

    @Test(enabled = true, dataProvider = "data")
    public void testDiffToString(DifferenceTest test) {
        logger.warn("Test tree1: " + (test.tree1 == null ? "null" : test.tree1.toOneLineString()));
        logger.warn("Test tree2: " + (test.tree2 == null ? "null" : test.tree2.toOneLineString()));
        logger.warn("Test expected diff : " + test.difference);
        Difference diff = new Difference(test.tree1, test.tree2);
        logger.warn("Observed diffs     : " + diff);
        Assert.assertEquals(diff.toString(), test.difference, "Observed diff string " + diff + " not equal to expected difference string " + test.difference);

    }
}