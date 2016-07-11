package org.broadinstitute.hellbender.utils.tsv;

import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;

/**
 * Unit tests for {@link TableColumnCollection}.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public class TableColumnCollectionUnitTest extends BaseTest {

    private static int[] COUNT_TEST = new int[]{1, 3, 10, 101};

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testCheckNamesOnANullArray() {
        TableColumnCollection.checkNames(null, GATKException::new);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testCheckNamesOnNullExceptionFactory() {
        TableColumnCollection.checkNames(new String[]{"col1"}, null);
    }

    @Test(expectedExceptions = GATKException.class)
    public void testCheckNamesOnNoColumnArray() {
        TableColumnCollection.checkNames(new String[0], GATKException::new);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testCheckNamesOnNullContainingArray() {
        TableColumnCollection.checkNames(new String[]{"col1", null, "col3"}, GATKException::new);
    }

    @Test(expectedExceptions = GATKException.class)
    public void testCheckNamesOnCommentPrefixedArray() {
        TableColumnCollection.checkNames(new String[]{"#col1", "col2", "col3"}, GATKException::new);
    }

    @Test(expectedExceptions = GATKException.class)
    public void testCheckNamesWithRepeats() {
        TableColumnCollection.checkNames(new String[]{"col1", "col2", "col1"}, GATKException::new);
    }

    @Test(dataProvider= "correctColumnNamesData")
    public void testCheckNamesOnCorrectData(final String[] names) {
        Assert.assertSame(TableColumnCollection.checkNames(names, IllegalArgumentException::new), names);
    }

    @Test(dataProvider = "correctColumnNamesData")
    public void testArrayCreation(final String[] names) {
        Assert.assertNotNull(new TableColumnCollection(names));
    }

    @Test(dataProvider = "correctColumnNamesData")
    public void testIterableCreation(final String[] names) {
        Assert.assertNotNull(new TableColumnCollection(Arrays.asList(names)));
    }

    @Test(dataProvider = "correctColumnNamesData", dependsOnMethods = "testArrayCreation")
    public void testColumnCount(final String[] names) {
        Assert.assertEquals(new TableColumnCollection(names).columnCount(), names.length);
    }

    @Test(dataProvider = "correctColumnNamesData", dependsOnMethods = "testArrayCreation")
    public void testMatches(final String[] names) {
        final String[] testNames = names.clone();
        final TableColumnCollection subject = new TableColumnCollection(testNames);
        for (int i = 0; i < testNames.length; i++) {
            Assert.assertTrue(subject.matches(i, names[i]));
            Assert.assertFalse(subject.matches(i + 1, names[i]));
            if (names.length > 1) {
                Assert.assertFalse(subject.matches(i, i == testNames.length - 1 ? names[i - 1] : names[i + 1]));
            }
            for (int j = i; j < testNames.length; j++) {
                final String[] subset = Arrays.copyOfRange(testNames, i, j + 1);
                Assert.assertTrue(subject.matchesAll(i, subset));
                Assert.assertFalse(subject.matchesAll(i + 1, subset));
            }
            Assert.assertTrue(subject.matchesAll(i));
        }
    }

    @Test(dataProvider = "correctColumnNamesData", dependsOnMethods = "testArrayCreation")
    public void testContains(final String[] names) {
        final String[] testNames = names.clone();
        final TableColumnCollection subject = new TableColumnCollection(testNames);
        for (int i = 0; i < testNames.length; i++) {
            Assert.assertTrue(subject.contains(names[i]));
            Assert.assertFalse(subject.contains(names[i] + "_no_there"));
            for (int j = i; j < testNames.length; j++) {
                final String[] subset = Arrays.copyOfRange(testNames, i, j + 1);
                Assert.assertTrue(subject.containsAll(subset));
                subset[subset.length >> 1] = subset[0] + "_no_there";
                Assert.assertFalse(subject.containsAll(subset));
            }
            Assert.assertTrue(subject.matchesAll(i));
        }
    }

    @Test(dataProvider = "correctColumnNamesData", dependsOnMethods = "testArrayCreation")
    public void testNameAt(final String[] names) {
        final String[] testNames = names.clone();
        final TableColumnCollection subject = new TableColumnCollection(testNames);
        for (int i = 0; i < testNames.length; i++) {
            Assert.assertEquals(subject.nameAt(i), names[i]);
        }
    }

    @Test(dataProvider = "correctColumnNamesData", dependsOnMethods = "testArrayCreation")
    public void testNames(final String[] names) {
        final String[] testNames = names.clone();
        final TableColumnCollection subject = new TableColumnCollection(testNames);
        Assert.assertEquals(subject.names(), Arrays.asList(names));
    }

    @Test(dataProvider = "correctColumnNamesData", dependsOnMethods = "testArrayCreation")
    public void testIndexOf(final String[] names) {
        final String[] testNames = names.clone();
        final TableColumnCollection subject = new TableColumnCollection(testNames);
        for (int i = 0; i < testNames.length; i++) {
            Assert.assertEquals(subject.indexOf(names[i]), i);
        }
    }

    @Test(dataProvider = "correctColumnNamesData", dependsOnMethods = "testArrayCreation")
    public void testMatchesExactly(final String[] names) {
        final String[] testNames = names.clone();
        final TableColumnCollection subject = new TableColumnCollection(testNames);
        Assert.assertTrue(subject.matchesExactly(testNames));

        // try changing a single name.
        for (int i = 0; i < testNames.length; i++) {
            final String[] namesChanged = Arrays.copyOf(testNames, testNames.length);
            namesChanged[i] = testNames[i] + "_changed";
            Assert.assertFalse(subject.matchesExactly(namesChanged));
        }

        // try switching two names.
        if (testNames.length > 1) {
            for (int i = 0; i < testNames.length; i++) {
                for (int j = i + 1; j < testNames.length; j++) {
                    final String[] namesChanged = Arrays.copyOf(testNames, testNames.length);
                    namesChanged[i] = testNames[j];
                    namesChanged[j] = testNames[i];
                    Assert.assertFalse(subject.matchesExactly(namesChanged));
                }
            }
        }

        // try removing one name from the end.
        final String[] shortenNames = Arrays.copyOf(testNames, testNames.length - 1);
        Assert.assertFalse(subject.matchesExactly(shortenNames));

        // try vs the empty string
        Assert.assertFalse(subject.matchesExactly());

        // try vs an extended name array.
        final String[] longerNames = Arrays.copyOf(testNames, testNames.length + 1);
        longerNames[testNames.length] = "extra_col";
        Assert.assertFalse(subject.matchesExactly(longerNames));

        testNames[0] = testNames[0] + "_changed";
        Assert.assertFalse(subject.matchesExactly(testNames));
    }

    @Test(dataProvider = "correctColumnNamesData", dependsOnMethods = "testArrayCreation")
    public void testContainsExactly(final String[] names) {
        final String[] testNames = names.clone();
        final TableColumnCollection subject = new TableColumnCollection(testNames);
        Assert.assertTrue(subject.containsExactly(testNames));

        // try changing a single name.
        for (int i = 0; i < testNames.length; i++) {
            final String[] namesChanged = Arrays.copyOf(testNames, testNames.length);
            namesChanged[i] = testNames[i] + "_changed";
            Assert.assertFalse(subject.containsExactly(namesChanged));
        }

        // try switching two names.
        if (testNames.length > 1) {
            for (int i = 0; i < testNames.length; i++) {
                for (int j = i + 1; j < testNames.length; j++) {
                    final String[] namesChanged = Arrays.copyOf(testNames, testNames.length);
                    namesChanged[i] = testNames[j];
                    namesChanged[j] = testNames[i];
                    Assert.assertTrue(subject.containsExactly(namesChanged));
                }
            }
        }

        // try removing one name from the end.
        final String[] shortenNames = Arrays.copyOf(testNames, testNames.length - 1);
        Assert.assertFalse(subject.containsExactly(shortenNames));

        // try vs the empty string
        Assert.assertFalse(subject.containsExactly());

        // try vs an extended name array.
        final String[] longerNames = Arrays.copyOf(testNames, testNames.length + 1);
        longerNames[testNames.length] = "extra_col";
        Assert.assertFalse(subject.containsExactly(longerNames));

        testNames[0] = testNames[0] + "_changed";
        Assert.assertFalse(subject.containsExactly(testNames));
    }


    @DataProvider(name = "correctColumnNamesData")
    public Object[][] correctColumnNamesData() {
        final Object[][] result = new Object[COUNT_TEST.length][];

        for (int i = 0; i < COUNT_TEST.length; i++) {
            final String[] names = new String[COUNT_TEST[i]];
            for (int j = 0; j < COUNT_TEST[i]; j++) {
                names[j] = "col_" + (j + 1);
            }
            result[i] = new Object[]{names};
        }
        return result;
    }
}
