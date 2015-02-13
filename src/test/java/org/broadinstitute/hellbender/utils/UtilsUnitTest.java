/*
* Copyright (c) 2012 The Broad Institute
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.hellbender.utils;

import org.apache.commons.io.FileUtils;
import org.apache.commons.lang3.ArrayUtils;
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
                final List<Integer> left = new LinkedList<Integer>();
                for ( int i = 0; i < leftSize; i++ ) left.add(i);
                final List<Integer> total = new LinkedList<Integer>();
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
    public void testDupStringNoChars() {
        String duped = Utils.dupString('a',0);
        Assert.assertEquals(duped.length(), 0, "dupString did not produce zero-length string");
    }

    @Test
    public void testDupStringOneChar() {
        String duped = Utils.dupString('b',1);
        Assert.assertEquals(duped.length(), 1, "dupString did not produce single character string");
        Assert.assertEquals(duped.charAt(0), 'b', "dupString character was incorrect");
    }

    @Test
    public void testXor() {
        Assert.assertEquals(Utils.xor(false, false), false, "xor F F failed");
        Assert.assertEquals(Utils.xor(false, true), true, "xor F T failed");
        Assert.assertEquals(Utils.xor(true, false), true, "xor T F failed");
        Assert.assertEquals(Utils.xor(true, true), false, "xor T T failed");
    }

    @Test
    public void testDupStringMultiChar() {
        String duped = Utils.dupString('c',5);
        Assert.assertEquals(duped.length(), 5, "dupString did not produce five character string");
        Assert.assertEquals(duped,"ccccc","dupString string was incorrect");
    }

    @Test
    public void testJoinMap() {
        Map<String,Integer> map = new LinkedHashMap<String,Integer>();
        map.put("one",1);
        map.put("two",2);
        String joined = Utils.joinMap("-",";",map);
        Assert.assertTrue("one-1;two-2".equals(joined));
    }

    @Test
    public void testJoinMapLargerSet() {
        Map<String,Integer> map = new LinkedHashMap<String,Integer>();
        map.put("one",1);
        map.put("two",2);
        map.put("three",1);
        map.put("four",2);
        map.put("five",1);
        map.put("six",2);
        String joined = Utils.joinMap("-",";",map);
        Assert.assertTrue("one-1;two-2;three-1;four-2;five-1;six-2".equals(joined));
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
    public void testListFromPrimitives(){
        int[] ints = {1,2,3,4};
        Assert.assertEquals(Utils.listFromPrimitives(ints), Arrays.asList(ArrayUtils.toObject(ints)));
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
    }

    @Test
    public void testConcat() {
        final String s1 = "A";
        final String s2 = "CC";
        final String s3 = "TTT";
        final String s4 = "GGGG";
        Assert.assertEquals(new String(Utils.concat()), "");
        Assert.assertEquals(new String(Utils.concat(s1.getBytes())), s1);
        Assert.assertEquals(new String(Utils.concat(s1.getBytes(), s2.getBytes())), s1 + s2);
        Assert.assertEquals(new String(Utils.concat(s1.getBytes(), s2.getBytes(), s3.getBytes())), s1 + s2 + s3);
        Assert.assertEquals(new String(Utils.concat(s1.getBytes(), s2.getBytes(), s3.getBytes(), s4.getBytes())), s1 + s2 + s3 + s4);
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
            final int[] valuesClone = values.clone();
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
            final double[] valuesClone = values.clone();
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

    @Test
    public void testCalcMD5() throws Exception {
        final File source = new File(publicTestDir + "exampleFASTA.fasta");
        final String sourceMD5 = "d0d4cbffece546c231fabd0103f54592";

        final byte[] sourceBytes = IOUtils.readFileIntoByteArray(source);
        Assert.assertEquals(Utils.calcMD5(sourceBytes), sourceMD5);

        final String sourceString = FileUtils.readFileToString(source);
        Assert.assertEquals(Utils.calcMD5(sourceString), sourceMD5);
    }

    @Test
    public void testLongestCommonOps() {
        for ( int prefixLen = 0; prefixLen < 20; prefixLen++ ) {
            for ( int extraSeq1Len = 0; extraSeq1Len < 10; extraSeq1Len++ ) {
                for ( int extraSeq2Len = 0; extraSeq2Len < 10; extraSeq2Len++ ) {
                    for ( int max = 0; max < 50; max++ ) {
                        final String prefix = Utils.dupString("A", prefixLen);
                        final int expected = Math.min(prefixLen, max);

                        {
                            final String seq1 = prefix + Utils.dupString("C", extraSeq1Len);
                            final String seq2 = prefix + Utils.dupString("G", extraSeq1Len);
                            Assert.assertEquals(Utils.longestCommonPrefix(seq1.getBytes(), seq2.getBytes(), max), expected, "LongestCommonPrefix failed: seq1 " + seq1 + " seq2 " + seq2 + " max " + max);
                        }

                        {
                            final String seq1 = Utils.dupString("C", extraSeq1Len) + prefix;
                            final String seq2 = Utils.dupString("G", extraSeq1Len) + prefix;
                            Assert.assertEquals(Utils.longestCommonSuffix(seq1.getBytes(), seq2.getBytes(), max), expected, "longestCommonSuffix failed: seq1 " + seq1 + " seq2 " + seq2 + " max " + max);
                        }
                    }
                }
            }
        }
    }

    @DataProvider(name = "trim")
    public Object[][] createTrimTestData() {
        List<Object[]> tests = new ArrayList<Object[]>();

        final String s = "AAAA";
        for ( int front = 0; front < s.length(); front++ ) {
            for ( int back = 0; back < s.length(); back++ ) {
                if ( front + back <= s.length() )
                    tests.add(new Object[]{s, front, back});
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "trim", enabled = true)
    public void testTrim(final String s, final int frontTrim, final int backTrim) {
        Assert.assertEquals(s.length() - frontTrim - backTrim, Utils.trimArray(s.getBytes(), frontTrim, backTrim).length);
    }

    @Test(dataProvider = "equalRangeData", enabled = true)
    public void testEqualRange(final byte[] array1, final byte[] array2, final int offset1, final int offset2, final int length, final boolean expected) {
        Assert.assertEquals(Utils.equalRange(array1,offset1,array2,offset2,length),expected);
        Assert.assertTrue(Utils.equalRange(array1,offset1,array1,offset1,length));
        Assert.assertTrue(Utils.equalRange(array2,offset2,array2,offset2,length));

    }

    @DataProvider(name = "equalRangeData")
    public Object[][] equalRangeData() {
        return new Object[][] {
                new Object[] { new byte[0] , new byte[0], 0, 0, 0, true},
                new Object[]  {      "ABCF".getBytes(), "BC".getBytes(), 1,0,2, true },
                new Object[]  { "ABCF".getBytes(), "".getBytes(), 1,0,0, true },
                new Object[]  { "ABCF".getBytes(), "ACBF".getBytes(), 0,0, 4, false}
        };

    }

    @Test(dataProvider = "skimArrayData")
    public void testSkimArray(final String original, final String remove) {
        final StringBuilder resultBuilder = new StringBuilder();
        final boolean[] removeBoolean = new boolean[remove.length()];
        for (int i = 0; i < original.length(); i++)
            if (remove.charAt(i) == '1') {
                resultBuilder.append(original.charAt(i));
                removeBoolean[i] = false;
            } else
                removeBoolean[i] = true;

        final String expected = resultBuilder.toString();
        final byte[] resultBytes = Utils.skimArray(original.getBytes(),removeBoolean);
        final String resultString = new String(resultBytes);
        Assert.assertEquals(resultString,expected);
    }

    @DataProvider(name = "skimArrayData")
    public Object[][] skimArrayData() {
        return new Object[][] {
                {"romeo+juliette" , "11111111111111" },
                {"romeo+juliette" , "11111011111111" },
                {"romeo+juliette" , "00000011111111" },
                {"romeo+juliette" , "11111100000000" },
                {"romeo+juliette" , "11111011111111" },
                {"romeo+juliette" , "01111010000001" },
                {"romeo+juliette" , "01100110000110" },
                {"romeo+juliette" , "10101010101010" },
                {"romeo+juliette" , "01010101010101" },
                {"romeo+juliette" , "01111010111001" },
        };
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
}
