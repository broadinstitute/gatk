package org.broadinstitute.hellbender.utils;


import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableMap;
import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.google.common.primitives.Ints;
import htsjdk.samtools.util.Log.LogLevel;
import org.apache.commons.io.FileUtils;
import org.apache.logging.log4j.Level;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import static java.util.Arrays.asList;
import static org.testng.Assert.assertEquals;

/**
 * Testing framework for general purpose utilities class.
 *
 */
public final class UtilsUnitTest extends GATKBaseTest {

    @Test
    public void testForceJVMLocaleToUSEnglish() {

        // Set locale to Canada
        Locale.setDefault(Locale.CANADA);

        // Force Locale to US English
        Utils.forceJVMLocaleToUSEnglish();

        // Get the current locale
        Locale l = Locale.getDefault();

        Assert.assertEquals(l, Locale.US);
    }

    @Test
    public void testConcat() {
        check(ImmutableList.of(ImmutableList.of()), ImmutableList.of());
        check(ImmutableList.of(ImmutableList.of("a")), ImmutableList.of("a"));
        check(ImmutableList.of(ImmutableList.of("a", "b")), ImmutableList.of("a", "b"));
        check(ImmutableList.of(ImmutableList.of("a"), ImmutableList.of("b")), ImmutableList.of("a", "b"));
        check(ImmutableList.of(ImmutableList.of(), ImmutableList.of("a"), ImmutableList.of("b")), ImmutableList.of("a", "b"));
        check(ImmutableList.of(ImmutableList.of("a"), ImmutableList.of(), ImmutableList.of("b")), ImmutableList.of("a", "b"));
        check(ImmutableList.of(ImmutableList.of("a"), ImmutableList.of("b"), ImmutableList.of()), ImmutableList.of("a", "b"));
        check(ImmutableList.of(ImmutableList.of("a", "b"), ImmutableList.of("c", "d")),
                ImmutableList.of("a", "b", "c", "d"));
    }

    private <T> void check(List<? extends Iterable<T>> input, List<T> expected) {
        assertEquals(Lists.newArrayList(Utils.concatIterators(input.iterator())), expected);
    }

    @Test
    public void testTransformParallel() {
        final Iterator<Integer> integers = Utils.transformParallel(ImmutableList.of(5, 4, 3, 2, 1).iterator(), i -> {
            try { Thread.sleep(i * 100); } catch (InterruptedException e) { }
            return i;
        }, 2);
        assertEquals(Lists.newArrayList(integers), ImmutableList.of(5, 4, 3, 2, 1));
    }

    @Test
    public void testIteratorConcat() throws Exception {
        final List<Integer> ints1 = Arrays.asList(0, 1, 2);
        final List<Integer> ints2 = Arrays.asList(3, 4, 5);
        final Iterator<Integer> it = Utils.concatIterators(Arrays.asList(ints1, ints2).iterator());
        final List<Integer> lst = Lists.newArrayList(it);
        Assert.assertEquals(lst, Arrays.asList(0,1,2, 3, 4, 5));
    }

    @Test
    public void testTransform() throws Exception {
        final Iterator<String> it= Arrays.asList("1", "2", "3").iterator();
        final Iterator<Integer> objectIterator = Utils.transformParallel(it, Integer::parseInt, 1);
        final List<Integer> lst = Lists.newArrayList(objectIterator);
        Assert.assertEquals(lst, Arrays.asList(1,2,3));
    }

    @Test
    public void testTransformInParallel() throws Exception {
        final Iterator<String> it= Arrays.asList("1", "2", "3").iterator();
        final Iterator<Integer> objectIterator = Utils.transformParallel(it, Integer::parseInt, 2);
        final List<Integer> lst = Lists.newArrayList(objectIterator);
        Assert.assertEquals(lst, Arrays.asList(1,2,3));
    }

    @Test
    public void testXor()  {
        Assert.assertEquals(Utils.xor(false, false), false);
        Assert.assertEquals(Utils.xor(false, true),  true);
        Assert.assertEquals(Utils.xor(true, false),  true);
        Assert.assertEquals(Utils.xor(true, true),   false);
    }

    @Test
    public void testRepeatBytes() throws Exception {
        Assert.assertEquals(Utils.repeatBytes((byte)112, 4), new byte[]{112,112,112,112});
    }

    @Test
    public void testRepeatChars() throws Exception {
        Assert.assertEquals(Utils.repeatChars('a', 4), new byte[]{'a','a','a','a'});
    }

    @Test
    public void testMakePermutations(){
//        * if objects = [A, B, C]
//        * if N = 1 => [[A], [B], [C]]
//        * if N = 2 => [[A, A], [B, A], [C, A], [A, B], [B, B], [C, B], [A, C], [B, C], [C, C]]

        final List<List<String>> p1 = Utils.makePermutations(asList("A", "B", "C"), 1, false);
        final List<List<String>> r1 = Utils.makePermutations(asList("A", "B", "C"), 1, true);

        final List<List<String>> p2 = Utils.makePermutations(asList("A", "B", "C"), 2, false);
        final List<List<String>> r2 = Utils.makePermutations(asList("A", "B", "C"), 2, true);

        Assert.assertEquals(asList(asList("A"), asList("B"), asList("C")), p1);
        Assert.assertEquals(asList(asList("A"), asList("B"), asList("C")), r1);

        Assert.assertEquals(asList(asList("B", "A"), asList("C", "A"), asList("A", "B"), asList("C", "B"), asList("A", "C"), asList("B", "C")), p2);
        Assert.assertEquals(asList(asList("A", "A"), asList("B", "A"), asList("C", "A"), asList("A", "B"), asList("B", "B"), asList("C", "B"), asList("A", "C"), asList("B", "C"), asList("C", "C")), r2);

    }

    @Test
    public void testJoinCollection() {
        Assert.assertEquals("0-1-2", Utils.join("-", asList(0, 1, 2)));
        Assert.assertEquals("0", Utils.join("-", asList(0)));
        Assert.assertEquals("", Utils.join("-", asList()));
    }

    @Test
    public void testCons() {
        final List<Integer> list = asList(0, 1, 2, 3);
        Assert.assertEquals(asList(6, 0, 1, 2, 3), Utils.cons(6, list));
    }

    @Test
    public void testDupBytes() {
        Assert.assertEquals(new byte[]{(byte)'a',(byte)'a',(byte)'a'}, Utils.dupBytes((byte) 'a', 3));
    }

    @Test
    public void testDupChar() {
        Assert.assertEquals("aaa", Utils.dupChar('a', 3));
        Assert.assertEquals("     ", Utils.dupChar(' ', 5));
    }

    @Test
    public void testAppend() {
        for ( int leftSize : asList(0, 1, 2, 3) ) {
            for ( final int rightSize : asList(0, 1, 2) ) {
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
        Assert.assertEquals(Utils.warnUserLines(message), new ArrayList<>(asList(
                "**********************************************************************",
                "* WARNING:",
                "* ",
                "* Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do",
                "* eiusmod tempor incididunt ut labore et dolore magna aliqua.",
                "**********************************************************************")));
    }


    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testJoinNullInts() {
        int[] nullints = null;
        Utils.join(",", nullints);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testJoinNullDoubles() {
        double[] emptydbl = null;
        Utils.join(",", emptydbl);
    }

    @Test
    public void testJoin(){
        int[] ints = {1,2,3,4};
        Assert.assertEquals(Utils.join(",", ints),"1,2,3,4");

        double[] dbls = {1.0,2.0,3.0,4.0};
        Assert.assertEquals(Utils.join(",",dbls), "1.0,2.0,3.0,4.0");

        Assert.assertEquals(Utils.join(",", new Object[] {}),"");
        Assert.assertEquals(Utils.join(",", new Object[] { true , -12, "Blah", this.getClass() }),
                "true,-12,Blah," + this.getClass().toString());
        Assert.assertEquals(Utils.join(",", true, -13, "Blah", this.getClass()),"true,-13,Blah," + this.getClass().toString());
        Assert.assertEquals(Utils.join(",", Boolean.TRUE), "true");
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

        final String sourceString = FileUtils.readFileToString(source, StandardCharsets.UTF_8);
        Assert.assertEquals(Utils.calcMD5(sourceString), sourceMD5);

        Assert.assertEquals(Utils.calculateFileMD5(source), sourceMD5);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testNonNullThrows(){
        final Object o = null;
        Utils.nonNull(o);
    }

    @Test
    public void testNonNullDoesNotThrow(){
        final Object o = new Object();
        Assert.assertSame(Utils.nonNull(o), o);
    }

    @Test(expectedExceptions = IllegalArgumentException.class, expectedExceptionsMessageRegExp = "^The exception message$")
    public void testNonNullWithMessageThrows() {
        Utils.nonNull(null, "The exception message");
    }

    @Test
    public void testNonNullWithMessageReturn() {
        final Object testObject = new Object();
        Assert.assertSame(Utils.nonNull(testObject, "some message"), testObject);
    }

    @Test
    public void testNonEmpty(){
        final Collection<String> notEmpty = Collections.singletonList("string1");
        Assert.assertSame(Utils.nonEmpty(notEmpty), notEmpty);
    }

    @DataProvider(name= "emptyAndNull")
    public Object[][] emptyAndNull(){
        return new Object[][] {
                {Collections.emptyList()},
                {null}
        };
    }

    @Test(expectedExceptions = IllegalArgumentException.class, dataProvider = "emptyAndNull")
    public void testNonEmptyThrows(Collection<?> collection) {
        Utils.nonEmpty(collection);
    }

    @Test(expectedExceptions = IllegalArgumentException.class, dataProvider = "emptyAndNull")
    public void testNonEmptyThrowsWithMessage(Collection<?> collection) {
        Utils.nonEmpty(collection, "some message");
    }


    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testCheckForDuplicatesWithDuplicates() {
        final List<String> strings = Arrays.asList("A", "B", "B");
        Utils.checkForDuplicatesAndReturnSet(strings, "Uh-oh");
    }

    @Test
    public void testCheckForDuplicatesWithoutDuplicates() {
        final List<String> strings = Arrays.asList("A", "B", "C", "E", "D");
        final Set<String> set = Utils.checkForDuplicatesAndReturnSet(strings, "Uh-oh");
        Assert.assertEquals(strings.stream().sorted().collect(Collectors.toList()), set.stream().sorted().collect(Collectors.toList()));
    }

    @Test(dataProvider = "successfulValidIndexData")
    public void testSuccessfulValidIndex(final int index, final int length) {
        final int actualIndex = Utils.validIndex(index, length);
        Assert.assertSame(actualIndex, index);
    }

    @Test(dataProvider = "unsuccessfulValidIndexData", expectedExceptions = IllegalArgumentException.class)
    public void testUnsuccessfulValidIndex(final int index, final int length) {
        Utils.validIndex(index, length);
    }

    @DataProvider(name = "successfulValidIndexData")
    public Object[][] successfulValidIndexData() {
        return new Object[][]{
                {0, 10}, {1, 2}
        };
    }

    @DataProvider(name = "unsuccessfulValidIndexData")
    public Object[][] unsuccessfulValidIndexData() {
        return new Object[][]{
                {0, 0}, {10, 10}, {11, 5}, {-1, 10}, {0, -10}
        };
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

    @DataProvider(name = "booleanOccurrencesData")
    public Object[][] booleanOccurrencesTestData() {
        return new Object[][]{
                new Object[]{ true, new boolean[] { }, 0 },
                new Object[]{ false, new boolean[] { }, 0 },
                new Object[]{ true, new boolean[] { true }, 1 },
                new Object[]{ false, new boolean[] { true }, 0 },
                new Object[]{ true, new boolean[] { false }, 0 },
                new Object[]{ true, new boolean[] { true, true, true}, 3},
                new Object[]{ false, new boolean[] { true, true, true}, 0},
                new Object[]{ true, new boolean[] { true, true, false}, 2},
        };
    }

    @Test(dataProvider = "booleanOccurrencesData")
    public void testCountBooleanOccurrences(boolean element, boolean[] array, int expected) {
        Assert.assertEquals(Utils.countBooleanOccurrences(element, array), expected);
    }

    /**
     * Test setting the global logging level for Picard, Log4j, MinLog and java.util.logging.
     *
     * Note that there are three very similar, but not identical, logging level enums from different namespaces
     * being used here. The one used by Picard (and Hellbender --verbosity) of type "Log.LogLevel", the parallel
     * one used by log4j of type "Level", and the one used by java.utils.logging.
     *
     */
    @Test
    public void testSetLoggingLevel() {
        // Query and cache the initial  Log4j level in place at the start of the tests so we can restore it at the
        // end of the tests. Also, since we're QUERYING the Log4j level, but we're SETTING the level using the
        // LoggingUtils API, we also need to verify here that the initial level is one of the narrower set of levels
        // that is supported by LoggingUtils, since those are the only ones we can restore through the LoggingUtils API.
        Level initialLevel = logger.getLevel();
        boolean goodInitialLevel =
                initialLevel == Level.DEBUG ||
                initialLevel == Level.WARN ||
                initialLevel == Level.ERROR ||
                initialLevel == Level.INFO;
        Assert.assertTrue(goodInitialLevel);

        // Set and test each supported logging level in turn
        LoggingUtils.setLoggingLevel(LogLevel.DEBUG);
        Assert.assertTrue(logger.getLevel() == Level.DEBUG);

        LoggingUtils.setLoggingLevel(LogLevel.WARNING);
        Assert.assertTrue(logger.getLevel() == Level.WARN);

        LoggingUtils.setLoggingLevel(LogLevel.ERROR);
        Assert.assertTrue(logger.getLevel() == Level.ERROR);

        LoggingUtils.setLoggingLevel(LogLevel.INFO);
        Assert.assertTrue(logger.getLevel() == Level.INFO);

        // Restore the logging level back to the original level in place at the beginning of the test
        LoggingUtils.setLoggingLevel(LoggingUtils.levelFromLog4jLevel(initialLevel));
        Assert.assertTrue(logger.getLevel() == initialLevel);
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
        Assert.assertEquals(resultString, expected);
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

    @DataProvider(name = "testEqualRange")
    public Object[][] testEqualRange(){
        final Object[][] result = {
            {"ATGATGATGATG".getBytes(), 0, 3, "ATG".getBytes(), true},
            {"ATGATGATGATG".getBytes(), 3, 6, "ATG".getBytes(), true},
            {"ATGATGATGATG".getBytes(), 9, 12, "ATG".getBytes(), true},

            {"ATGATGATGATG".getBytes(), 0, 1, "G".getBytes(), false},
            {"ATGATGATCATG".getBytes(), 0, 2, "AT".getBytes(), true},
            {"ATGATGATCATG".getBytes(), 2, 4, "AT".getBytes(), false},

            {"T".getBytes(), 0, 1, "T".getBytes(), true},

            {"CCCCCCCC".getBytes(), 0, 3, "CCC".getBytes(), true},
            {"CCCCCCCC".getBytes(), 5, 8, "CCC".getBytes(), true},

            {"GATAT".getBytes(), 3, 5, "AT".getBytes(), true},
            {"GATAT".getBytes(), 1, 3, "AT".getBytes(), true},

            {"ATGATGATGATG".getBytes(), 11, 12, "G".getBytes(), true},
            {"ATGATGATGATG".getBytes(), 10, 11, "G".getBytes(), false},
            {"ATGATGATCATG".getBytes(), 10, 12, "AT".getBytes(), false},
        };
        return result;
    }

    @Test(dataProvider = "testEqualRange")
    public void testEqualRange(byte[] arr1, int from, int to, byte[] arr2, boolean expected) throws Exception {
        //Note 'from' is inclusive, 'to' is exclusive
        Assert.assertEquals(Utils.equalRange(arr1, from, arr2, 0, to - from), expected);
    }

    @Test
    public void testLastIndexOfQueryTooLong() {
        final String reference = "AAAA";
        final String query     = "AAAAAAA";

        final int result = Utils.lastIndexOf(reference.getBytes(), query.getBytes());
        final int expected = reference.lastIndexOf(query);
        Assert.assertEquals(result, expected);
    }

    @Test
    public void testLastIndexOfLastBoundaries() {
        final String reference = "AAAACCCCTTTTGGGG";

        // match right boundary of reference
        String query = "TGGGG";
        int result = Utils.lastIndexOf(reference.getBytes(), query.getBytes());
        int expected = reference.lastIndexOf(query);
        Assert.assertEquals(result, expected);

        // match left boundary of reference
        query = "AAAAC";
        result = Utils.lastIndexOf(reference.getBytes(), query.getBytes());
        expected = reference.lastIndexOf(query);
        Assert.assertEquals(result, expected);
    }

    private void randomByteString(Random rng, byte[] bytes) {
        for (int i = 0; i < bytes.length; i++) {
            bytes[i] = (byte)(rng.nextInt(94) + 32);
        }
    }

    @Test
    public void testLastIndexOfRandom() {
        final int num_tests = 100;
        final int referenceLength = 1000;
        final int queryLength = 100;
        
        byte [] reference = new byte[referenceLength];
        byte [] query = new byte[queryLength];

        final Random rng = Utils.getRandomGenerator();
        
        for (int i = 0; i < num_tests; i++) {
            randomByteString(rng, reference);
            randomByteString(rng, query);

            // add query to reference at a random location for 75% of the tests
            if (i % 4 > 0) {
                final int index = rng.nextInt(referenceLength - queryLength);
                for (int j = 0; j < queryLength; j++) {
                    reference[index+j] = query[j];
                }
            }
            
            final int result = Utils.lastIndexOf(reference, query);
            final int expected = new String(reference).lastIndexOf(new String(query));
            Assert.assertEquals(result, expected);
        }
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testListFromPrimitivesNull() throws Exception {
        Utils.listFromPrimitives(null);
    }

    @Test
    public void testListFromPrimitivesEmpty() throws Exception {
        Assert.assertTrue(Utils.listFromPrimitives(new int[0]).isEmpty());
    }

    @Test
    public void testListFromPrimitivesNoneEmpty() throws Exception {
        Assert.assertEquals(Utils.listFromPrimitives(new int[]{1,2}), Arrays.asList(1,2));
    }

    @DataProvider
    public Object[][] getNonNullCollections(){
        final List<String> someValues = Arrays.asList("some", "values");
        return new Object[][]{
                {Collections.emptyList()},
                {Collections.emptySet()},
                {someValues},
                {new HashSet<>(someValues)},
                {new TreeSet<>(someValues)},
        };
    }

    @DataProvider
    public Object[][] getCollectionsWithNulls(){
        return new Object[][]{
                {null},
                {Arrays.asList("something", null)},
                {Sets.newHashSet("something", null)},
        };
    }

    @Test(dataProvider = "getNonNullCollections")
    public void testContainsNoNull(Collection<?> collection){
        Utils.containsNoNull(collection, "bad");
    }

    @Test(dataProvider = "getCollectionsWithNulls", expectedExceptions = IllegalArgumentException.class)
    public void testContainsNull( Collection<?> collection){
        Utils.containsNoNull(collection, "This was expected");
    }

    @Test
    public void testNonEmptyString() {
        Assert.assertThrows(IllegalArgumentException.class, () -> Utils.nonEmpty((String)null, "deliberately empty"));
        Assert.assertThrows(IllegalArgumentException.class, () -> Utils.nonEmpty("", "deliberately empty" ));
        Utils.nonEmpty("this is not empty");
    }

    @Test
    public void testStreamFromIterable() {
        final int[] array = new int[] {1,3,5,7,9};
        Assert.assertTrue(Arrays.equals(array, Utils.stream(Ints.asList(array)).mapToInt(n -> n).toArray()));
    }

    @Test
    public void testStreamFromIterator() {
        final int[] array = new int[] {1,3,5,7,9};
        Assert.assertEquals(array, Utils.stream(Ints.asList(array).iterator()).mapToInt(n -> n).toArray());
    }

    @DataProvider
    public Object[][] getCollectionsForDuplicatedItemsTest() {
        return new Object[][] {
                {Collections.emptyList(), Collections.emptySet()},
                {Collections.singletonList("U"), Collections.emptySet()},
                {Arrays.asList("U1", "U2"), Collections.emptySet()},
                {Arrays.asList("D", "D"), Collections.singleton("D")},
                {Arrays.asList("U", "D", "D"), Collections.singleton("D")},
                {Arrays.asList("D", "U", "D"), Collections.singleton("D")},
                {Arrays.asList("D", "D", "U"), Collections.singleton("D")},
                {Arrays.asList("D1", "D1", "D2", "D2"), new HashSet<>(Arrays.asList("D1", "D2"))},
                {Arrays.asList("D1", "D2", "D1", "D2"), new HashSet<>(Arrays.asList("D1", "D2"))},
                {Arrays.asList("D1", "D2", "D2", "D1"), new HashSet<>(Arrays.asList("D1", "D2"))}
        };
    }

    @Test(dataProvider = "getCollectionsForDuplicatedItemsTest")
    public void testGetDuplicatedItems(final Collection<?> collection, final Set<?> duplicated) {
        final Set<?> result = Utils.getDuplicatedItems(collection);
        Assert.assertEquals(result, duplicated);
    }

    @DataProvider
    public Iterator<Object[]> provideDataForTestUtilsSplitString() {

        final String stringData = "The quick fox jumped over the lazy brown dog.  " +
                "Arma virumque cano, Troiae qui primus ab oris " +
                "Italiam, fato profugus, Laviniaque venit " +
                "litora, multum ille et terris iactatus et alto " +
                "vi superum saevae memorem Iunonis ob iram; " +
                "multa quoque et bello passus, testConfigurationSorting conderet urbem, " +
                "inferretque deos Latio, genus unde Latinum, " +
                "Albanique patres, atque altae moenia Romae.  " +
                "Musa, mihi causas memora, quo numine laeso, " +
                "quidve dolens, regina deum tot volvere casus " +
                "insignem pietate virum, tot adire labores " +
                "impulerit. Tantaene animis caelestibus irae9";

        final List<String> wordsToSplitOn = Arrays.stream(stringData.split(" "))
                .map(ssss -> (ssss.contains("?")) ? ssss.replace("?", "\\?") : ssss)
                .collect(Collectors.toList());

        final List<String> repeatedSubstringsToSplitOn = Arrays.asList( "us", "is", "it", "et", " et", " et ", "que ", " mult" );

        final String allDelimiterTestString = "::::::::::::::::::::::::::::::";
        final String frontDelimiterTestString = "::::::::::1:23::::::::::456:7890";
        final String middleDelimiterTestString = "1:23::::::::::456:7890";
        final String backDelimiterTestString = "1:23::::::::::456:7890::::::::::";
        final String fullDelimiterTestString = "::::::::::1:23::::::::::456:7890::::::::::";

        final List<Object[]> testCases = new ArrayList<>();
        testCases.addAll(
            Arrays.asList(
                new Object[] { "", "" },
                new Object[] { ":", "" },
                new Object[] { "SOME MORE TESTS", "" },
                new Object[] { ":", ":" },
                new Object[] { stringData, ""  },
                new Object[] { stringData, "1" },
                new Object[] { stringData, "a" },
                new Object[] { stringData, "b" },
                new Object[] { stringData, "c" },
                new Object[] { stringData, "d" },
                new Object[] { stringData, "e" },
                new Object[] { stringData, "f" },
                new Object[] { stringData, "g" },
                new Object[] { stringData, "h" },
                new Object[] { stringData, "i" },
                new Object[] { stringData, "j" },
                new Object[] { stringData, "k" },
                new Object[] { stringData, "l" },
                new Object[] { stringData, "m" },
                new Object[] { stringData, "n" },
                new Object[] { stringData, "o" },
                new Object[] { stringData, "p" },
                new Object[] { stringData, "q" },
                new Object[] { stringData, "r" },
                new Object[] { stringData, "s" },
                new Object[] { stringData, "t" },
                new Object[] { stringData, "u" },
                new Object[] { stringData, "v" },
                new Object[] { stringData, "w" },
                new Object[] { stringData, "x" },
                new Object[] { stringData, "y" },
                new Object[] { stringData, "z" },
                new Object[] { stringData, " " },
                new Object[] { stringData, "T" },
                new Object[] { stringData, "9" },
                new Object[] { allDelimiterTestString,    ":" },
                new Object[] { frontDelimiterTestString,  ":" },
                new Object[] { middleDelimiterTestString, ":" },
                new Object[] { backDelimiterTestString ,  ":" },
                new Object[] { fullDelimiterTestString ,  ":" },
                new Object[] { allDelimiterTestString.replace(":", "TS"),    "TS" },
                new Object[] { frontDelimiterTestString.replace(":", "TS"),  "TS" },
                new Object[] { middleDelimiterTestString.replace(":", "TS"), "TS" },
                new Object[] { backDelimiterTestString.replace(":", "TS") ,  "TS" },
                new Object[] { fullDelimiterTestString.replace(":", "TS") ,  "TS" }
            )
        );

        // Create test cases for words:
        for ( final String delim : wordsToSplitOn ) {
            testCases.add( new Object[] { stringData, delim } );
        }

        // Create test cases for repeated substrings:
        for ( final String delim : repeatedSubstringsToSplitOn ) {
            testCases.add( new Object[] { stringData, delim } );
        }

        return testCases.iterator();
    }

    @DataProvider
    private Iterator<Object[]> provideDataForTestUtilsSplitStringExhaustively() {

        // Length that we want to test through:
        final int maxStringLength = 7;

        // Create single character delimiter strings:
        final String singleCharDelimiter = "o";
        final List<String> singleCharTestStrings = Arrays.asList("X", singleCharDelimiter);
        final List<List<String>> exhaustiveListsForSingleChar = Utils.makePermutations( singleCharTestStrings, maxStringLength, true );

        // Create multi-character delimiter strings:
        final String multiCharDelimiter = "oz";

        // Must include the individual characters of the multi-char delimiter here
        // so we can pass them into Utils.makePermutations to create the permutations of strings to split.
        final List<String> multiCharTestStrings = Arrays.asList("X", "o", "z", multiCharDelimiter);
        final List<List<String>> exhaustiveListsForMultiChar = Utils.makePermutations( multiCharTestStrings, maxStringLength, true );

        final List<Object[]> testCases = new ArrayList<>();

        // Add single-char cases:
        for ( final List<String> testCase : exhaustiveListsForSingleChar ) {
            testCases.add( new Object[] { String.join( "", testCase ), singleCharDelimiter } );
        }

        // Add multi-char cases:
        for ( final List<String> testCase : exhaustiveListsForMultiChar ) {
            testCases.add( new Object[] { String.join( "", testCase ), multiCharDelimiter } );
        }

        return testCases.iterator();
    }

    private void exhaustiveStringSplitHelper(final String str,
                                             final String delimiter) {
        List<String> splitStrings;
        if ( delimiter.length() == 1 ) {
            splitStrings = Utils.split( str, delimiter.charAt(0) );
            Assert.assertEquals( splitStrings, Arrays.asList(str.split(delimiter)) );
        }

        splitStrings = Utils.split(str, delimiter);
        Assert.assertEquals( splitStrings, Arrays.asList(str.split(delimiter)) );
    }

    @Test(dataProvider = "provideDataForTestUtilsSplitString")
    public void testUtilsSplitString( final String str, final String delimiter ) {
        exhaustiveStringSplitHelper(str, delimiter);
    }

    @Test(dataProvider = "provideDataForTestUtilsSplitStringExhaustively")
    public void testUtilsSplitStringExhaustively( final String str, final String delimiter ) {
        exhaustiveStringSplitHelper(str, delimiter);
    }

    @DataProvider
    public Object[][] provideGetReverseValueToListMap() {
        return new Object[][]{
                {ImmutableMap.of("Foo", Arrays.asList(1, 2, 3)),
                        ImmutableMap.of(1, Sets.newHashSet("Foo"), 2, Sets.newHashSet("Foo"), 3, Sets.newHashSet("Foo"))},
                {ImmutableMap.of("Foo", Arrays.asList(1, 2, 3), "Baz", Arrays.asList(1, 2)),
                        ImmutableMap.of(1, Sets.newHashSet("Foo", "Baz"), 2, Sets.newHashSet("Foo", "Baz"), 3, Sets.newHashSet("Foo"))},

                // Let's throw in some nulls
                {Collections.unmodifiableMap(Stream.of(
                    new AbstractMap.SimpleEntry<>(null, Arrays.asList("one", "two", "three")),
                    new AbstractMap.SimpleEntry<>("two", Arrays.asList("one", "two")))
                    .collect(Collectors.toMap(AbstractMap.SimpleEntry::getKey, AbstractMap.SimpleEntry::getValue))),
                 Collections.unmodifiableMap(Stream.of(
                    new AbstractMap.SimpleEntry<>("one", Sets.newHashSet(null, "two")),
                    new AbstractMap.SimpleEntry<>("two", Sets.newHashSet(null, "two")),
                    new AbstractMap.SimpleEntry<>("three", Sets.newHashSet((Object) null)))
                    .collect(Collectors.toMap(AbstractMap.SimpleEntry::getKey, AbstractMap.SimpleEntry::getValue)))},

                // Let's throw in some nulls again
                {Collections.unmodifiableMap(Stream.of(
                    new AbstractMap.SimpleEntry<>(null, Arrays.asList(null, "two", "three")),
                    new AbstractMap.SimpleEntry<>("two", Arrays.asList("one", "two")))
                    .collect(Collectors.toMap(AbstractMap.SimpleEntry::getKey, AbstractMap.SimpleEntry::getValue))),
                 Collections.unmodifiableMap(Stream.of(
                    new AbstractMap.SimpleEntry<>("one", Sets.newHashSet( "two")),
                    new AbstractMap.SimpleEntry<>("two", Sets.newHashSet(null, "two")),
                    new AbstractMap.SimpleEntry<>("three", Sets.newHashSet((Object) null)),
                    new AbstractMap.SimpleEntry<>(null, Sets.newHashSet((Object) null)))
                    .collect(Collectors.toMap(AbstractMap.SimpleEntry::getKey, AbstractMap.SimpleEntry::getValue)))},

                    // Let's throw in some nulls again and a non-string
                {Collections.unmodifiableMap(Stream.of(
                    new AbstractMap.SimpleEntry<>(null, Arrays.asList(null, 2, "three")),
                    new AbstractMap.SimpleEntry<>("two", Arrays.asList("one", "two")))
                    .collect(Collectors.toMap(AbstractMap.SimpleEntry::getKey, AbstractMap.SimpleEntry::getValue))),
                 Collections.unmodifiableMap(Stream.of(
                    new AbstractMap.SimpleEntry<>("one", Sets.newHashSet( "two")),
                    new AbstractMap.SimpleEntry<>("two", Sets.newHashSet( "two")),
                    new AbstractMap.SimpleEntry<>("three", Sets.newHashSet((Object) null)),
                    new AbstractMap.SimpleEntry<>(null, Sets.newHashSet((Object) null)),
                    new AbstractMap.SimpleEntry<>(2, Sets.newHashSet((Object) null)))
                    .collect(Collectors.toMap(AbstractMap.SimpleEntry::getKey, AbstractMap.SimpleEntry::getValue)))}
        };
    }

    @Test(dataProvider = "provideGetReverseValueToListMap")
    public <T,U> void testGetReverseValueToListMap(final Map<T, List<U>> input,  final Map<U, Set<T>> gtOutput) {
        Assert.assertEquals(Utils.getReverseValueToListMap(input), gtOutput);
    }

    @DataProvider
    public Iterator<Object[]> provideDataForTestUtilsFormatting() {
        final List<Object[]> testCases = new ArrayList<>();
        testCases.addAll(Arrays.asList(
                new Object[]{1, 2, "50.00", "0.50"},
                new Object[]{0, 3, "0.00", "0.00"},
                new Object[]{1, 3, "33.33", "0.33"},
                new Object[]{50, 3000, "1.67", "0.02"},
                new Object[]{50, 0, "NA", "NA"},
                new Object[]{0, 0, "NA", "NA"}
        ));

        return testCases.iterator();
    }

    @Test(dataProvider = "provideDataForTestUtilsFormatting")
    public void testFormattedPctAndRatio(final long input1, final long input2, final String formattedPct, final String formattedRatio) {
        Assert.assertEquals(Utils.formattedPercent(input1, input2), formattedPct);
        Assert.assertEquals(Utils.formattedRatio(input1, input2), formattedRatio);
    }

    @DataProvider(name="provideDataForTestFilterCollectionByExpressions")
    public Object[][] provideDataForTestFilterCollectionByExpressions() {
        return new Object[][] {
                new Object[] { new LinkedHashSet<>(Arrays.asList("a", "ab", "abc")), Arrays.asList("a"), true, new LinkedHashSet<>(Arrays.asList("a")) },
                new Object[] { new LinkedHashSet<>(Arrays.asList("a", "ab", "abc")), Arrays.asList("a"), false, new LinkedHashSet<>(Arrays.asList("a", "ab", "abc")) },
                new Object[] { new LinkedHashSet<>(Arrays.asList("a", "ab", "abc")), Arrays.asList("b"), true, Collections.EMPTY_SET },
                new Object[] { new LinkedHashSet<>(Arrays.asList("a", "ab", "abc")), Arrays.asList("b"), false, new LinkedHashSet<>(Arrays.asList("ab", "abc")) },
                new Object[] { new LinkedHashSet<>(Arrays.asList("a", "ab", "abc")), Arrays.asList("a", "b"), true, new LinkedHashSet<>(Arrays.asList("a")) },
                new Object[] { new LinkedHashSet<>(Arrays.asList("a", "ab", "abc")), Arrays.asList("a", "b"), false, new LinkedHashSet<>(Arrays.asList("a", "ab", "abc")) },
                new Object[] { new LinkedHashSet<>(Arrays.asList("a", "ab", "abc")), Arrays.asList("a", "ab"), true, new LinkedHashSet<>(Arrays.asList("a", "ab")) },
                new Object[] { new LinkedHashSet<>(Arrays.asList("a", "ab", "abc")), Arrays.asList("a", "ab"), false, new LinkedHashSet<>(Arrays.asList("a", "ab", "abc")) },
                new Object[] { new LinkedHashSet<>(Arrays.asList("a", "ab", "abc")), Arrays.asList(".*b.*"), true, Collections.EMPTY_SET },
                new Object[] { new LinkedHashSet<>(Arrays.asList("a", "ab", "abc")), Arrays.asList(".*b.*"), false, new LinkedHashSet<>(Arrays.asList("ab", "abc") )},
                new Object[] { new LinkedHashSet<>(Arrays.asList("a", "ab", "abc")), Arrays.asList(".*"), true, Collections.EMPTY_SET },
                new Object[] { new LinkedHashSet<>(Arrays.asList("a", "ab", "abc")), Arrays.asList(".*"), false, new LinkedHashSet<>(Arrays.asList("a", "ab", "abc") )}
        };
    }

    @Test(dataProvider = "provideDataForTestFilterCollectionByExpressions")
    public void testTestFilterCollectionByExpressions(Set<String> values, Collection<String> filters, boolean exactMatch, Set<String> expected) {
        Set<String> actual = Utils.filterCollectionByExpressions(values, filters, exactMatch);
        Assert.assertEquals(actual, expected);
    }
}
