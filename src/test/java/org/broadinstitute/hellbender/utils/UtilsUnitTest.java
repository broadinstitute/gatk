package org.broadinstitute.hellbender.utils;

import org.apache.commons.io.FileUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.*;

/**
 * Testing framework for general purpose utilities class.
 *
 */
public class UtilsUnitTest extends BaseTest {

    @Test
    public void testAppend() {
        for ( int leftSize : Arrays.asList(0, 1, 2, 3) ) {
            for ( final int rightSize : Arrays.asList(0, 1, 2) ) {
                final List<Integer> left = new LinkedList<>();
                for ( int i = 0; i < leftSize; i++ ) left.add(i);
                final List<Integer> total = new LinkedList<>();
                for ( int i = 0; i < leftSize + rightSize; i++ ) total.add(i);

                if ( rightSize == 0 )
                    Assert.assertEquals(Utils.append(left), total);
                if ( rightSize == 1 )
                    Assert.assertEquals(Utils.append(left, leftSize), total);
                if ( rightSize == 2 )
                    Assert.assertEquals(Utils.append(left, leftSize, leftSize + 1), total);
            }
        }

    }

    @Test
    public void testWarnUserLines(){
        String message = "Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut" +
                " labore et dolore magna aliqua.";
        Assert.assertEquals(Utils.warnUserLines(message), new ArrayList<>(Arrays.asList(
                "**********************************************************************",
                "* WARNING:",
                "* ",
                "* Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do",
                "* eiusmod tempor incididunt ut labore et dolore magna aliqua.",
                "**********************************************************************")));
    }

    @Test
    public void testJoin(){
        int[] ints = {1,2,3,4};
        int[] nullints = null;
        Assert.assertEquals(Utils.join(",", ints),"1,2,3,4");
        Assert.assertEquals(Utils.join(",", nullints), "");

        double[] dbls = {1.0,2.0,3.0,4.0};
        double[] emptydbl = null;
        Assert.assertEquals(Utils.join(",", emptydbl), "");
        Assert.assertEquals(Utils.join(",",dbls), "1.0,2.0,3.0,4.0");

        Assert.assertEquals(Utils.join(",", new Object[] {}),"");
        Assert.assertEquals(Utils.join(",", new Object[] { true , -12, "Blah", this.getClass() }),
                "true,-12,Blah," + this.getClass().toString());
        Assert.assertEquals(Utils.join(",", true, -13, "Blah", this.getClass()),"true,-13,Blah," + this.getClass().toString());
        Assert.assertEquals(Utils.join(",", Boolean.TRUE),"true");
    }

    @Test
    public void testEscapeExpressions() {
        String[] expected, actual;

        expected = new String[] {"one", "two", "three"};
        actual = Utils.escapeExpressions("one two three");
        Assert.assertEquals(actual, expected);
        actual = Utils.escapeExpressions(" one two three");
        Assert.assertEquals(actual, expected);
        actual = Utils.escapeExpressions("one two three ");
        Assert.assertEquals(actual, expected);
        actual = Utils.escapeExpressions(" one two three ");
        Assert.assertEquals(actual, expected);
        actual = Utils.escapeExpressions("  one  two  three  ");
        Assert.assertEquals(actual, expected);

        expected = new String[] {"one", "two", "three four", "five", "six"};
        actual = Utils.escapeExpressions("one two 'three four' five six");
        Assert.assertEquals(actual, expected);
        actual = Utils.escapeExpressions(" one two 'three four' five six");
        Assert.assertEquals(actual, expected);
        actual = Utils.escapeExpressions("one two 'three four' five six ");
        Assert.assertEquals(actual, expected);
        actual = Utils.escapeExpressions(" one two 'three four' five six ");
        Assert.assertEquals(actual, expected);
        actual = Utils.escapeExpressions("  one  two  'three four'  five  six  ");
        Assert.assertEquals(actual, expected);

        expected = new String[] {"one two", "three", "four"};
        actual = Utils.escapeExpressions("'one two' three four");
        Assert.assertEquals(actual, expected);
        actual = Utils.escapeExpressions(" 'one two' three four");
        Assert.assertEquals(actual, expected);
        actual = Utils.escapeExpressions("'one two' three four ");
        Assert.assertEquals(actual, expected);
        actual = Utils.escapeExpressions(" 'one two' three four ");
        Assert.assertEquals(actual, expected);
        actual = Utils.escapeExpressions("  'one two'  three  four  ");
        Assert.assertEquals(actual, expected);

        expected = new String[] {"one", "two", "three four"};
        actual = Utils.escapeExpressions("one two 'three four'");
        Assert.assertEquals(actual, expected);
        actual = Utils.escapeExpressions(" one two 'three four'");
        Assert.assertEquals(actual, expected);
        actual = Utils.escapeExpressions("one two 'three four' ");
        Assert.assertEquals(actual, expected);
        actual = Utils.escapeExpressions(" one two 'three four' ");
        Assert.assertEquals(actual, expected);
        actual = Utils.escapeExpressions("  one  two  'three four'  ");
        Assert.assertEquals(actual, expected);
    }

    @Test(dataProvider = "asIntegerListData")
    public void testAsIntegerList(final int[] values) {
        if (values == null) {
            try {
                Utils.asList((int[]) null);
                Assert.fail("Should have thrown an exception");
            } catch (final IllegalArgumentException ex) {
                // good.
            }
        } else {
            final Random rdn = Utils.getRandomGenerator();
            final int[] valuesClone = Arrays.copyOf(values, values.length);
            final List<Integer> list = Utils.asList(valuesClone);
            Assert.assertNotNull(list);
            Assert.assertEquals(list.size(),values.length);
            for (int i = 0; i < values.length; i++)
                Assert.assertEquals((int) list.get(i),values[i]);
            for (int i = 0; i < values.length; i++)
                valuesClone[rdn.nextInt(values.length)] = rdn.nextInt(1000);
            for (int i = 0; i < values.length; i++)
                Assert.assertEquals((int) list.get(i),valuesClone[i]);
        }
    }

    @DataProvider(name = "asIntegerListData")
    public Object[][] asIntegerListData() {
        return new Object[][] {
                { null },
                {new int[0]},
                {new int[]{1, 2, 3, 4, 5}},
                {new int[]{2}},
                {new int[]{3,4}}
        };
    }

    @Test(dataProvider = "asDoubleListData")
    public void testAsDoubleList(final double[] values) {
        if (values == null) {
            try {
                Utils.asList((int[]) null);
                Assert.fail("Should have thrown an exception");
            } catch (final IllegalArgumentException ex) {
                // good.
            }
        } else {
            final Random rdn = Utils.getRandomGenerator();
            final double[] valuesClone = Arrays.copyOf(values, values.length);
            final List<Double> list = Utils.asList(valuesClone);
            Assert.assertNotNull(list);
            Assert.assertEquals(list.size(),values.length);
            for (int i = 0; i < values.length; i++)
                Assert.assertEquals(list.get(i),values[i]);
            for (int i = 0; i < values.length; i++)
                valuesClone[rdn.nextInt(values.length)] = rdn.nextDouble() * 1000;
            for (int i = 0; i < values.length; i++)
                Assert.assertEquals(list.get(i),valuesClone[i]);
        }
    }


    @DataProvider(name = "asDoubleListData")
    public Object[][] asDoubleListData() {
        return new Object[][] {
                { null },
                {new double[0]},
                {new double[]{1, 2, 3, 4, 5}},
                {new double[]{2}},
                {new double[]{3,4}},
                {new double[]{Double.NaN, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY}}
        };
    }

    @Test
    public void testCalcMD5() throws Exception {
        final File source = new File(publicTestDir + "exampleFASTA.fasta");
        final String sourceMD5 = "d0d4cbffece546c231fabd0103f54592";

        final byte[] sourceBytes = IOUtils.readFileIntoByteArray(source);
        Assert.assertEquals(Utils.calcMD5(sourceBytes), sourceMD5);

        final String sourceString = FileUtils.readFileToString(source);
        Assert.assertEquals(Utils.calcMD5(sourceString), sourceMD5);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testNonNullThrows(){
        Object o = null;
        Utils.nonNull(o);
    }

    @Test
    public void testNonNullDoesNotThrow(){
        Object o = new Object();
        Assert.assertSame(Utils.nonNull(o), o);
    }
}
