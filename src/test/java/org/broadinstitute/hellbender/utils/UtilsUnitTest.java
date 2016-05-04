package org.broadinstitute.hellbender.utils;

import htsjdk.samtools.util.Log.LogLevel;
import org.apache.commons.io.FileUtils;
import org.apache.logging.log4j.Level;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.SkipException;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.*;

import static java.util.Arrays.asList;

/**
 * Testing framework for general purpose utilities class.
 *
 */
public final class UtilsUnitTest extends BaseTest {

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

        final String sourceString = FileUtils.readFileToString(source);
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


    @Test
    public void testSuccessfulRegularReadableFileCheck() {
        final File expectedFile = createTempFile("Utils-RRFC-test",".txt");
        final File actualFile = Utils.regularReadableUserFile(expectedFile);
        Assert.assertSame(actualFile, expectedFile);
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

    @Test(dataProvider = "unsuccessfulRegularReadableFileCheckData",
            expectedExceptions = UserException.CouldNotReadInputFile.class)
    public void testUnsuccessfulRegularReadableFileCheck(final File file) {
        if (file == null){
            throw new SkipException("cannot make a file unreadable (maybe you're running as root)");
        }
        Utils.regularReadableUserFile(file);
    }

    @DataProvider(name = "unsuccessfulRegularReadableFileCheckData")
    public Object[][] unsuccessfulRegularReadableFileCheckData() {
        final File directory = createTempDir("Utils-RRFCD-Dir");
        final File nonExistingFile = createTempFile("Utils-RRFCD-NoFile", ".file");
        nonExistingFile.delete();
        final File nonReadable = createTempFile("Utils-RRFCD-NoReadable", ".file");

        //Note: if this test suite is run as root (eg in a Docker image), setting readable to false may fail.
        // So we check it here and skip this test if we can't make a file non readable.
        final boolean successSetNotReadable = nonReadable.setReadable(false) && !nonReadable.canRead();
        if (successSetNotReadable) {
            return new Object[][]{{directory}, {nonExistingFile}, {nonReadable}};
        } else {
            //pass in a special value to be able to throw a SkipException and make the same number of tests
            logger.debug("cannot make a file unreadable (maybe you're running as root)");
            return new Object[][]{{directory}, {nonExistingFile}, {null}};
        }
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
     * Test setting the global logging level for Picard and Log4j and java.util.logging.
     *
     * Note that there are three very similar, but not identical, logging level enums from different namespaces
     * being used here. The one used by Picard (and Hellbender --verbosity) of type "Log.LogLevel", the parallel
     * one used by log4j of type "Level", and the one used by java.utils.logging.
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
}
