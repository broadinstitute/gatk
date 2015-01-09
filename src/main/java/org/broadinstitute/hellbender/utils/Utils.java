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

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.lang.reflect.Array;
import java.math.BigInteger;
import java.net.InetAddress;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.*;

public class Utils {
    private Utils(){}

    /**
     *  Static random number generator and seed.
     */
    private static final long GATK_RANDOM_SEED = 47382911L;
    private static Random randomGenerator = new Random(GATK_RANDOM_SEED);
    public static Random getRandomGenerator() { return randomGenerator; }
    public static void resetRandomGenerator() { randomGenerator.setSeed(GATK_RANDOM_SEED); }
    public static void resetRandomGenerator(long seed) { randomGenerator.setSeed(seed); }

    private static final int TEXT_WARNING_WIDTH = 68;
    private static final String TEXT_WARNING_PREFIX = "* ";
    private static final String TEXT_WARNING_BORDER = dupString('*', TEXT_WARNING_PREFIX.length() + TEXT_WARNING_WIDTH);
    private static final char ESCAPE_CHAR = '\u001B';
    // ASCII codes for making text blink
    public static final String TEXT_BLINK = ESCAPE_CHAR + "[5m";
    public static final String TEXT_RESET = ESCAPE_CHAR + "[m";

    /** our log, which we want to capture anything from this class */
    private static Logger logger = LogManager.getLogger(Utils.class);

    public static final float JAVA_DEFAULT_HASH_LOAD_FACTOR = 0.75f;

    /**
     * Boolean xor operation.  Only true if x != y.
     *
     * @param x a boolean
     * @param y a boolean
     * @return true if x != y
     */
    public static boolean xor(final boolean x, final boolean y) {
        return x != y;
    }

    /**
     * Calculates the optimum initial size for a hash table given the maximum number
     * of elements it will need to hold. The optimum size is the smallest size that
     * is guaranteed not to result in any rehash/table-resize operations.
     *
     * @param maxElements  The maximum number of elements you expect the hash table
     *                     will need to hold
     * @return             The optimum initial size for the table, given maxElements
     */
    public static int optimumHashSize ( int maxElements ) {
        return (int)(maxElements / JAVA_DEFAULT_HASH_LOAD_FACTOR) + 2;
    }

    public static <T> List<T> cons(final T elt, final List<T> l) {
        List<T> l2 = new ArrayList<T>();
        l2.add(elt);
        if (l != null) l2.addAll(l);
        return l2;
    }

    public static void warnUser(final String msg) {
        warnUser(logger, msg);
    }

    public static void warnUser(final Logger logger, final String msg) {
        for (final String line: warnUserLines(msg))
            logger.warn(line);
    }

    public static List<String> warnUserLines(final String msg) {
        List<String> results = new ArrayList<>();
        results.add(String.format(TEXT_WARNING_BORDER));
        results.add(String.format(TEXT_WARNING_PREFIX + "WARNING:"));
        results.add(String.format(TEXT_WARNING_PREFIX));
        prettyPrintWarningMessage(results, msg);
        results.add(String.format(TEXT_WARNING_BORDER));
        return results;
    }

    /**
     * pretty print the warning message supplied
     *
     * @param results the pretty printed message
     * @param message the message
     */
    private static void prettyPrintWarningMessage(final List<String> results, final String message) {
        for (final String line: message.split("\\r?\\n")) {
            final StringBuilder builder = new StringBuilder(line);
            while (builder.length() > TEXT_WARNING_WIDTH) {
                int space = getLastSpace(builder, TEXT_WARNING_WIDTH);
                if (space <= 0) space = TEXT_WARNING_WIDTH;
                results.add(String.format("%s%s", TEXT_WARNING_PREFIX, builder.substring(0, space)));
                builder.delete(0, space + 1);
            }
            results.add(String.format("%s%s", TEXT_WARNING_PREFIX, builder));
        }
    }

    /**
     * Returns the last whitespace location in string, before width characters.
     * @param message The message to break.
     * @param width The width of the line.
     * @return The last whitespace location.
     */
    private static int getLastSpace(final CharSequence message, int width) {
        final int length = message.length();
        int stopPos = width;
        int currPos = 0;
        int lastSpace = -1;
        boolean inEscape = false;
        while (currPos < stopPos && currPos < length) {
            final char c = message.charAt(currPos);
            if (c == ESCAPE_CHAR) {
                stopPos++;
                inEscape = true;
            } else if (inEscape) {
                stopPos++;
                if (Character.isLetter(c))
                    inEscape = false;
            } else if (Character.isWhitespace(c)) {
                lastSpace = currPos;
            }
            currPos++;
        }
        return lastSpace;
    }

    /**
     * join the key value pairs of a map into one string, i.e. myMap = [A->1,B->2,C->3] with a call of:
     * joinMap("-","*",myMap) -> returns A-1*B-2*C-3
     *
     * Be forewarned, if you're not using a map that is aware of the ordering (i.e. HashMap instead of LinkedHashMap)
     * the ordering of the string you get back might not be what you expect! (i.e. C-3*A-1*B-2 vrs A-1*B-2*C-3)
     *
     * @param keyValueSeperator the string to seperate the key-value pairs
     * @param recordSeperator the string to use to seperate each key-value pair from other key-value pairs
     * @param map the map to draw from
     * @param <L> the map's key type
     * @param <R> the map's value type
     * @return a string representing the joined map
     */
    public static <L,R> String joinMap(String keyValueSeperator, String recordSeperator, Map<L,R> map) {
        if (map.size() < 1) { return null; }
        String joinedKeyValues[] = new String[map.size()];
        int index = 0;
        for (L key : map.keySet()) {
            joinedKeyValues[index++] = String.format("%s%s%s",key.toString(),keyValueSeperator,map.get(key).toString());
        }
        return join(recordSeperator,joinedKeyValues);
    }

    /**
     * Splits a String using indexOf instead of regex to speed things up.
     *
     * @param str the string to split.
     * @param delimiter the delimiter used to split the string.
     * @return an array of tokens.
     */
    public static ArrayList<String> split(String str, String delimiter) {
        return split(str, delimiter, 10);
    }

    /**
     * Splits a String using indexOf instead of regex to speed things up.
     *
     * @param str the string to split.
     * @param delimiter the delimiter used to split the string.
     * @param expectedNumTokens The number of tokens expected. This is used to initialize the ArrayList.
     * @return an array of tokens.
     */
    public static ArrayList<String> split(String str, String delimiter, int expectedNumTokens) {
        final ArrayList<String> result =  new ArrayList<String>(expectedNumTokens);

        int delimiterIdx = -1;
        do {
            final int tokenStartIdx = delimiterIdx + 1;
            delimiterIdx = str.indexOf(delimiter, tokenStartIdx);
            final String token = (delimiterIdx != -1 ? str.substring(tokenStartIdx, delimiterIdx) : str.substring(tokenStartIdx) );
            result.add(token);
        } while( delimiterIdx != -1 );

        return result;
    }


    /**
     * join an array of strings given a seperator
     * @param separator the string to insert between each array element
     * @param strings the array of strings
     * @return a string, which is the joining of all array values with the separator
     */
    public static String join(String separator, String[] strings) {
        return join(separator, strings, 0, strings.length);
    }

    public static String join(String separator, String[] strings, int start, int end) {
        if ((end - start) == 0) {
            return "";
        }
        StringBuilder ret = new StringBuilder(strings[start]);
        for (int i = start + 1; i < end; ++i) {
            ret.append(separator);
            ret.append(strings[i]);
        }
        return ret.toString();
    }

    public static String join(String separator, int[] ints) {
        if ( ints == null || ints.length == 0)
            return "";
        else {
            StringBuilder ret = new StringBuilder();
            ret.append(ints[0]);
            for (int i = 1; i < ints.length; ++i) {
                ret.append(separator);
                ret.append(ints[i]);
            }
            return ret.toString();
        }
    }

    /**
     * Returns a string of the values in joined by separator, such as A,B,C
     *
     * @param separator separator character
     * @param doubles   the array with values
     * @return a string with the values separated by the separator
     */
    public static String join(String separator, double[] doubles) {
        if ( doubles == null || doubles.length == 0)
            return "";
        else {
            StringBuilder ret = new StringBuilder();
            ret.append(doubles[0]);
            for (int i = 1; i < doubles.length; ++i) {
                ret.append(separator);
                ret.append(doubles[i]);
            }
            return ret.toString();
        }
    }

    /**
     * Returns a string of the form elt1.toString() [sep elt2.toString() ... sep elt.toString()] for a collection of
     * elti objects (note there's no actual space between sep and the elti elements).  Returns
     * "" if collection is empty.  If collection contains just elt, then returns elt.toString()
     *
     * @param separator the string to use to separate objects
     * @param objects a collection of objects.  the element order is defined by the iterator over objects
     * @param <T> the type of the objects
     * @return a non-null string
     */
    public static <T> String join(final String separator, final Collection<T> objects) {
        if (objects.isEmpty()) { // fast path for empty collection
            return "";
        } else {
            final Iterator<T> iter = objects.iterator();
            final T first = iter.next();

            if ( ! iter.hasNext() ) // fast path for singleton collections
                return first.toString();
            else { // full path for 2+ collection that actually need a join
                final StringBuilder ret = new StringBuilder(first.toString());
                while(iter.hasNext()) {
                    ret.append(separator);
                    ret.append(iter.next().toString());
                }
                return ret.toString();
            }
        }
    }

    /**
     * Returns a {@link List List&lt;Integer&gt;} representation of an primitive int array.
     * @param values the primitive int array to represent.
     * @return never code {@code null}. The returned list will be unmodifiable yet it will reflect changes in values in the original array yet
     *   you cannot change the values
     */
    public static List<Integer> asList(final int ... values) {
        if (values == null)
            throw new IllegalArgumentException("the input array cannot be null");
        return new AbstractList<Integer>() {

            @Override
            public Integer get(final int index) {
                return values[index];
            }

            @Override
            public int size() {
                return values.length;
            }
        };
    }

    /**
     * Returns a {@link List List&lt;Double&gt;} representation of an primitive double array.
     * @param values the primitive int array to represent.
     * @return never code {@code null}. The returned list will be unmodifiable yet it will reflect changes in values in the original array yet
     *   you cannot change the values.
     */
    public static List<Double> asList(final double ... values) {
        if (values == null)
            throw new IllegalArgumentException("the input array cannot be null");
        return new AbstractList<Double>() {

            @Override
            public Double get(final int index) {
                return values[index];
            }

            @Override
            public int size() {
                return values.length;
            }
        };
    }

    /**
     * Create a new list that contains the elements of left along with elements elts
     * @param left a non-null list of elements
     * @param elts a varargs vector for elts to append in order to left
     * @return A newly allocated linked list containing left followed by elts
     */
    public static <T> List<T> append(final List<T> left, T ... elts) {
        final List<T> l = new LinkedList<T>(left);
        l.addAll(Arrays.asList(elts));
        return l;
    }

    /**
     * Create a new string thats a n duplicate copies of s
     * @param s the string to duplicate
     * @param nCopies how many copies?
     * @return a string
     */
    public static String dupString(final String s, int nCopies) {
        if ( s == null || s.equals("") ) throw new IllegalArgumentException("Bad s " + s);
        if ( nCopies < 0 ) throw new IllegalArgumentException("nCopies must be >= 0 but got " + nCopies);

        final StringBuilder b = new StringBuilder();
        for ( int i = 0; i < nCopies; i++ )
            b.append(s);
        return b.toString();
    }

    public static String dupString(char c, int nCopies) {
        char[] chars = new char[nCopies];
        Arrays.fill(chars, c);
        return new String(chars);
    }

    public static byte[] dupBytes(byte b, int nCopies) {
        byte[] bytes = new byte[nCopies];
        Arrays.fill(bytes, b);
        return bytes;
    }

    // trim a string for the given character (i.e. not just whitespace)
    public static String trim(String str, char ch) {
        char[] array = str.toCharArray();


        int start = 0;
        while ( start < array.length && array[start] == ch )
            start++;

        int end = array.length - 1;
        while ( end > start && array[end] == ch )
            end--;

        return str.substring(start, end+1);
    }

    /**
     * Splits expressions in command args by spaces and returns the array of expressions.
     * Expressions may use single or double quotes to group any individual expression, but not both.
     * @param args Arguments to parse.
     * @return Parsed expressions.
     */
    public static String[] escapeExpressions(String args) {
        // special case for ' and " so we can allow expressions
        if (args.indexOf('\'') != -1)
            return escapeExpressions(args, "'");
        else if (args.indexOf('\"') != -1)
            return escapeExpressions(args, "\"");
        else
            return args.trim().split(" +");
    }

    /**
     * Splits expressions in command args by spaces and the supplied delimiter and returns the array of expressions.
     * @param args Arguments to parse.
     * @param delimiter Delimiter for grouping expressions.
     * @return Parsed expressions.
     */
    private static String[] escapeExpressions(String args, String delimiter) {
        String[] command = {};
        String[] split = args.split(delimiter);
        String arg;
        for (int i = 0; i < split.length - 1; i += 2) {
            arg = split[i].trim();
            if (arg.length() > 0) // if the unescaped arg has a size
                command = Utils.concatArrays(command, arg.split(" +"));
            command = Utils.concatArrays(command, new String[]{split[i + 1]});
        }
        arg = split[split.length - 1].trim();
        if (split.length % 2 == 1) // if the command ends with a delimiter
            if (arg.length() > 0) // if the last unescaped arg has a size
                command = Utils.concatArrays(command, arg.split(" +"));
        return command;
    }

    /**
     * Concatenates two String arrays.
     * @param A First array.
     * @param B Second array.
     * @return Concatenation of A then B.
     */
    public static String[] concatArrays(String[] A, String[] B) {
        String[] C = new String[A.length + B.length];
        System.arraycopy(A, 0, C, 0, A.length);
        System.arraycopy(B, 0, C, A.length, B.length);
        return C;
    }

    /**
     * Concatenates byte arrays
     * @return a concat of all bytes in allBytes in order
     */
    public static byte[] concat(final byte[] ... allBytes) {
        int size = 0;
        for ( final byte[] bytes : allBytes ) size += bytes.length;

        final byte[] c = new byte[size];
        int offset = 0;
        for ( final byte[] bytes : allBytes ) {
            System.arraycopy(bytes, 0, c, offset, bytes.length);
            offset += bytes.length;
        }

        return c;
    }

    /**
     * Appends String(s) B to array A.
     * @param A First array.
     * @param B Strings to append.
     * @return A with B(s) appended.
     */
    public static String[] appendArray(String[] A, String... B) {
        return concatArrays(A, B);
    }

    public static <T extends Comparable<T>> List<T> sorted(Collection<T> c) {
        return sorted(c, false);
    }

    public static <T extends Comparable<T>> List<T> sorted(Collection<T> c, boolean reverse) {
        List<T> l = new ArrayList<T>(c);
        Collections.sort(l);
        if ( reverse ) Collections.reverse(l);
        return l;
    }

    public static <T extends Comparable<T>, V> List<V> sorted(Map<T,V> c) {
        return sorted(c, false);
    }

    public static <T extends Comparable<T>, V> List<V> sorted(Map<T,V> c, boolean reverse) {
        List<T> t = new ArrayList<T>(c.keySet());
        Collections.sort(t);
        if ( reverse ) Collections.reverse(t);

        List<V> l = new ArrayList<V>();
        for ( T k : t ) {
            l.add(c.get(k));
        }
        return l;
    }

    /**
     * Reverse a byte array of bases
     *
     * @param bases  the byte array of bases
     * @return the reverse of the base byte array
     */
    static public byte[] reverse(byte[] bases) {
        byte[] rcbases = new byte[bases.length];

        for (int i = 0; i < bases.length; i++) {
            rcbases[i] = bases[bases.length - i - 1];
        }

        return rcbases;
    }

    static public <T> List<T> reverse(final List<T> l) {
        final List<T> newL = new ArrayList<T>(l);
        Collections.reverse(newL);
        return newL;
    }

    /**
     * Reverse an int array of bases
     *
     * @param bases  the int array of bases
     * @return the reverse of the base int array
     */
    static public int[] reverse(int[] bases) {
        int[] rcbases = new int[bases.length];

        for (int i = 0; i < bases.length; i++) {
            rcbases[i] = bases[bases.length - i - 1];
        }

        return rcbases;
    }

    /**
     * Reverse (NOT reverse-complement!!) a string
     *
     * @param bases  input string
     * @return the reversed string
     */
    static public String reverse(String bases) {
        return new String( reverse( bases.getBytes() )) ;
    }

    public static boolean isFlagSet(int value, int flag) {
        return ((value & flag) == flag);
    }

    /**
     * Helper utility that calls into the InetAddress system to resolve the hostname.  If this fails,
     * unresolvable gets returned instead.
     */
    public static String resolveHostname() {
        try {
            return InetAddress.getLocalHost().getCanonicalHostName();
        }
        catch (java.net.UnknownHostException uhe) { // [beware typo in code sample -dmw]
            return "unresolvable";
            // handle exception
        }
    }


    public static byte [] arrayFromArrayWithLength(byte[] array, int length) {
        byte [] output = new byte[length];
        for (int j = 0; j < length; j++)
            output[j] = array[(j % array.length)];
        return output;
    }

    public static void fillArrayWithByte(byte[] array, byte value) {
        for (int i=0; i<array.length; i++)
            array[i] = value;
    }


    /**
     * Returns the number of combinations represented by this collection
     * of collection of options.
     *
     * For example, if this is [[A, B], [C, D], [E, F, G]] returns 2 * 2 * 3 = 12
     */
    public static <T> int nCombinations(final Collection<T>[] options) {
        int nStates = 1;
        for ( Collection<T> states : options ) {
            nStates *= states.size();
        }
        return nStates;
    }

    public static <T> int nCombinations(final List<List<T>> options) {
        if ( options.isEmpty() )
            return 0;
        else {
            int nStates = 1;
            for ( Collection<T> states : options ) {
                nStates *= states.size();
            }
            return nStates;
        }
    }

    /**
     * Make all combinations of N size of objects
     *
     * if objects = [A, B, C]
     * if N = 1 => [[A], [B], [C]]
     * if N = 2 => [[A, A], [B, A], [C, A], [A, B], [B, B], [C, B], [A, C], [B, C], [C, C]]
     *
     * @param objects         list of objects
     * @param n               size of each combination
     * @param withReplacement if false, the resulting permutations will only contain unique objects from objects
     * @return a list with all combinations with size n of objects.
     */
    public static <T> List<List<T>> makePermutations(final List<T> objects, final int n, final boolean withReplacement) {
        final List<List<T>> combinations = new ArrayList<List<T>>();

        if ( n == 1 ) {
            for ( final T o : objects )
                combinations.add(Collections.singletonList(o));
        } else if (n > 1) {
            final List<List<T>> sub = makePermutations(objects, n - 1, withReplacement);
            for ( List<T> subI : sub ) {
                for ( final T a : objects ) {
                    if ( withReplacement || ! subI.contains(a) )
                        combinations.add(Utils.cons(a, subI));
                }
            }
        }

        return combinations;
    }

    /**
     * Convenience function that formats the novelty rate as a %.2f string
     *
     * @param known number of variants from all that are known
     * @param all number of all variants
     * @return a String novelty rate, or NA if all == 0
     */
    public static String formattedNoveltyRate(final int known, final int all) {
        return formattedPercent(all - known, all);
    }

    /**
     * Convenience function that formats the novelty rate as a %.2f string
     *
     * @param x number of objects part of total that meet some criteria
     * @param total count of all objects, including x
     * @return a String percent rate, or NA if total == 0
     */
    public static String formattedPercent(final long x, final long total) {
        return total == 0 ? "NA" : String.format("%.2f", (100.0*x) / total);
    }

    /**
     * Convenience function that formats a ratio as a %.2f string
     *
     * @param num  number of observations in the numerator
     * @param denom number of observations in the denumerator
     * @return a String formatted ratio, or NA if all == 0
     */
    public static String formattedRatio(final long num, final long denom) {
        return denom == 0 ? "NA" : String.format("%.2f", num / (1.0 * denom));
    }

    /**
     * Create a constant map that maps each value in values to itself
     */
    public static <T> Map<T, T> makeIdentityFunctionMap(Collection<T> values) {
        Map<T,T> map = new HashMap<T, T>(values.size());
        for ( final T value : values )
            map.put(value, value);
        return Collections.unmodifiableMap(map);
    }

    /**
     * Divides the input list into a list of sublists, which contains group size elements (except potentially the last one)
     *
     * list = [A, B, C, D, E]
     * groupSize = 2
     * result = [[A, B], [C, D], [E]]
     *
     */
    public static <T> List<List<T>> groupList(final List<T> list, final int groupSize) {
        if ( groupSize < 1 ) throw new IllegalArgumentException("groupSize >= 1");

        final List<List<T>> subLists = new LinkedList<List<T>>();
        int n = list.size();
        for ( int i = 0; i < n; i += groupSize ) {
            subLists.add(list.subList(i, Math.min(i + groupSize, n)));
        }
        return subLists;
    }

    /**
     * @see #calcMD5(byte[])
     */
    public static String calcMD5(final String s) {
        return calcMD5(s.getBytes());
    }

    /**
     * Calculate the md5 for bytes, and return the result as a 32 character string
     *
     * @param bytes the bytes to calculate the md5 of
     * @return the md5 of bytes, as a 32-character long string
     */
    public static String calcMD5(final byte[] bytes) {
        if ( bytes == null ) throw new IllegalArgumentException("bytes cannot be null");
        try {
            final byte[] thedigest = MessageDigest.getInstance("MD5").digest(bytes);
            final BigInteger bigInt = new BigInteger(1, thedigest);

            String md5String = bigInt.toString(16);
            while (md5String.length() < 32) md5String = "0" + md5String; // pad to length 32
            return md5String;
        }
        catch ( NoSuchAlgorithmException e ) {
            throw new IllegalStateException("MD5 digest algorithm not present");
        }
    }

    /**
     * Does big end with the exact sequence of bytes in suffix?
     *
     * @param big a non-null byte[] to test if it a prefix + suffix
     * @param suffix a non-null byte[] to test if it's a suffix of big
     * @return true if big is proper byte[] composed of some prefix + suffix
     */
    public static boolean endsWith(final byte[] big, final byte[] suffix) {
        if ( big == null ) throw new IllegalArgumentException("big cannot be null");
        if ( suffix == null ) throw new IllegalArgumentException("suffix cannot be null");
        return new String(big).endsWith(new String(suffix));
    }

    /**
     * Get the length of the longest common prefix of seq1 and seq2
     * @param seq1 non-null byte array
     * @param seq2 non-null byte array
     * @param maxLength the maximum allowed length to return
     * @return the length of the longest common prefix of seq1 and seq2, >= 0
     */
    public static int longestCommonPrefix(final byte[] seq1, final byte[] seq2, final int maxLength) {
        if ( seq1 == null ) throw new IllegalArgumentException("seq1 is null");
        if ( seq2 == null ) throw new IllegalArgumentException("seq2 is null");
        if ( maxLength < 0 ) throw new IllegalArgumentException("maxLength < 0 " + maxLength);

        final int end = Math.min(seq1.length, Math.min(seq2.length, maxLength));
        for ( int i = 0; i < end; i++ ) {
            if ( seq1[i] != seq2[i] )
                return i;
        }
        return end;
    }

    /**
     * Get the length of the longest common suffix of seq1 and seq2
     * @param seq1 non-null byte array
     * @param seq2 non-null byte array
     * @param maxLength the maximum allowed length to return
     * @return the length of the longest common suffix of seq1 and seq2, >= 0
     */
    public static int longestCommonSuffix(final byte[] seq1, final byte[] seq2, final int maxLength) {
        if ( seq1 == null ) throw new IllegalArgumentException("seq1 is null");
        if ( seq2 == null ) throw new IllegalArgumentException("seq2 is null");
        if ( maxLength < 0 ) throw new IllegalArgumentException("maxLength < 0 " + maxLength);

        final int end = Math.min(seq1.length, Math.min(seq2.length, maxLength));
        for ( int i = 0; i < end; i++ ) {
            if ( seq1[seq1.length - i - 1] != seq2[seq2.length - i - 1] )
                return i;
        }
        return end;
    }

    /**
     * Trim any number of bases from the front and/or back of an array
     *
     * @param seq                the sequence to trim
     * @param trimFromFront      how much to trim from the front
     * @param trimFromBack       how much to trim from the back
     * @return a non-null array; can be the original array (i.e. not a copy)
     */
    public static byte[] trimArray(final byte[] seq, final int trimFromFront, final int trimFromBack) {
        if ( trimFromFront + trimFromBack > seq.length )
            throw new IllegalArgumentException("trimming total is larger than the original array");

        // don't perform array copies if we need to copy everything anyways
        return  ( trimFromFront == 0 && trimFromBack == 0 ) ? seq : Arrays.copyOfRange(seq, trimFromFront, seq.length - trimFromBack);
    }

    /**
     * Simple wrapper for sticking elements of a int[] array into a List<Integer>
     * @param ar - the array whose elements should be listified
     * @return - a List<Integer> where each element has the same value as the corresponding index in @ar
     */
    public static List<Integer> listFromPrimitives(final int[] ar) {
        final ArrayList<Integer> lst = new ArrayList<>(ar.length);
        for ( final int d : ar ) {
            lst.add(d);
        }

        return lst;
    }

    /**
     * Compares sections from to byte arrays to verify whether they contain the same values.
     *
     * @param left first array to compare.
     * @param leftOffset first position of the first array to compare.
     * @param right second array to compare.
     * @param rightOffset first position of the second array to compare.
     * @param length number of positions to compare.
     *
     * @throws IllegalArgumentException if <ul>
     *     <li>either {@code left} or {@code right} is {@code null} or</li>
     *     <li>any off the offset or length combine point outside any of the two arrays</li>
     * </ul>
     * @return {@code true} iff {@code length} is 0 or all the bytes in both ranges are the same two-by-two.
     */
    public static boolean equalRange(final byte[] left, final int leftOffset, byte[] right, final int rightOffset, final int length) {
        if (left == null) throw new IllegalArgumentException("left cannot be null");
        if (right == null) throw new IllegalArgumentException("right cannot be null");
        if (length < 0) throw new IllegalArgumentException("the length cannot be negative");
        if (leftOffset < 0) throw new IllegalArgumentException("left offset cannot be negative");
        if (leftOffset + length > left.length) throw new IllegalArgumentException("length goes beyond end of left array");
        if (rightOffset < 0) throw new IllegalArgumentException("right offset cannot be negative");
        if (rightOffset + length > right.length) throw new IllegalArgumentException("length goes beyond end of right array");

        for (int i = 0; i < length; i++)
            if (left[leftOffset + i] != right[rightOffset + i])
                return false;
        return true;
    }

    /**
     * Skims out positions of an array returning a shorter one with the remaning positions in the same order.
     * @param original the original array to splice.
     * @param remove for each position in {@code original} indicates whether it should be spliced away ({@code true}),
     *               or retained ({@code false})
     *
     * @param <T> the array type.
     *
     * @throws IllegalArgumentException if either {@code original} or {@code remove} is {@code null},
     *    or {@code remove length is different to {@code original}'s}, or {@code original} is not in
     *    fact an array.
     *
     * @return never {@code null}.
     */
    public static <T> T skimArray(final T original, final boolean[] remove) {
        return skimArray(original,0,null,0,remove,0);
    }

    /**
     * Skims out positions of an array returning a shorter one with the remaning positions in the same order.
     *
     * <p>
     *     If the {@code dest} array provide is not long enough a new one will be created and returned with the
     *     same component type. All elements before {@code destOffset} will be copied from the input to the
     *     result array. If {@code dest} is {@code null}, a brand-new array large enough will be created where
     *     the position preceding {@code destOffset} will be left with the default value. The component type
     *     Will match the one of the {@code source} array.
     * </p>
     *
     * @param source the original array to splice.
     * @param sourceOffset the first position to skim.
     * @param dest the destination array.
     * @param destOffset the first position where to copy the skimed array values.
     * @param remove for each position in {@code original} indicates whether it should be spliced away ({@code true}),
     *               or retained ({@code false})
     * @param removeOffset the first position in the remove index array to consider.
     *
     * @param <T> the array type.
     *
     * @throws IllegalArgumentException if either {@code original} or {@code remove} is {@code null},
     *    or {@code remove length is different to {@code original}'s}, or {@code original} is not in
     *    fact an array.
     *
     * @return never {@code null}.
     */
    public static <T> T skimArray(final T source, final int sourceOffset, final T dest, final int destOffset, final boolean[] remove, final int removeOffset) {
        if (source == null)
            throw new IllegalArgumentException("the source array cannot be null");
        @SuppressWarnings("unchecked")
        final Class<T> sourceClazz = (Class<T>) source.getClass();

        if (!sourceClazz.isArray())
            throw new IllegalArgumentException("the source array is not in fact an array instance");
        final int length = Array.getLength(source) - sourceOffset;
        if (length < 0)
            throw new IllegalArgumentException("the source offset goes beyond the source array length");
        return skimArray(source,sourceOffset,dest,destOffset,remove,removeOffset,length);
    }

    /**
     * Skims out positions of an array returning a shorter one with the remaning positions in the same order.
     *
     * <p>
     *     If the {@code dest} array provide is not long enough a new one will be created and returned with the
     *     same component type. All elements before {@code destOffset} will be copied from the input to the
     *     result array. If {@code dest} is {@code null}, a brand-new array large enough will be created where
     *     the position preceding {@code destOffset} will be left with the default value. The component type
     *     Will match the one of the {@code source} array.
     * </p>
     *
     * @param source the original array to splice.
     * @param sourceOffset the first position to skim.
     * @param dest the destination array.
     * @param destOffset the first position where to copy the skimed array values.
     * @param remove for each position in {@code original} indicates whether it should be spliced away ({@code true}),
     *               or retained ({@code false})
     * @param removeOffset the first position in the remove index array to consider.
     * @param length the total number of position in {@code source} to consider. Thus only the {@code sourceOffset} to
     *               {@code sourceOffset + length - 1} region will be skimmed.
     *
     * @param <T> the array type.
     *
     * @throws IllegalArgumentException if either {@code original} or {@code remove} is {@code null},
     *    or {@code remove length is different to {@code original}'s}, or {@code original} is not in
     *    fact an array.
     *
     * @return never {@code null}.
     */
    public static <T> T skimArray(final T source, final int sourceOffset, final T dest, final int destOffset,
                                  final boolean[] remove, final int removeOffset, final int length) {
        if (source == null)
            throw new IllegalArgumentException("the source array cannot be null");
        if (remove == null)
            throw new IllegalArgumentException("the remove array cannot be null");
        if (sourceOffset < 0)
            throw new IllegalArgumentException("the source array offset cannot be negative");
        if (destOffset < 0)
            throw new IllegalArgumentException("the destination array offset cannot be negative");
        if (removeOffset < 0)
            throw new IllegalArgumentException("the remove array offset cannot be negative");
        if (length < 0)
            throw new IllegalArgumentException("the length provided cannot be negative");

        final int removeLength = Math.min(remove.length - removeOffset,length);

        if (removeLength < 0)
            throw new IllegalArgumentException("the remove offset provided falls beyond the remove array end");


        @SuppressWarnings("unchecked")
        final Class<T> sourceClazz = (Class<T>) source.getClass();

        if (!sourceClazz.isArray())
            throw new IllegalArgumentException("the source array is not in fact an array instance");

        final Class<T> destClazz = skimArrayDetermineDestArrayClass(dest, sourceClazz);

        final int sourceLength = Array.getLength(source);

        if (sourceLength < length + sourceOffset)
            throw new IllegalArgumentException("the source array is too small considering length and offset");

        // count how many positions are to be removed.

        int removeCount = 0;

        final int removeEnd = removeLength + removeOffset;
        for (int i = removeOffset; i < removeEnd; i++)
            if  (remove[i]) removeCount++;


        final int newLength = length - removeCount;


        @SuppressWarnings("unchecked")
        final T result = skimArrayBuildResultArray(dest, destOffset, destClazz, newLength);
        // No removals, just copy the whole thing.

        if (removeCount == 0)
            System.arraycopy(source,sourceOffset,result,destOffset,length);
        else if (length > 0) {  // if length == 0 nothing to do.
            int nextOriginalIndex = 0;
            int nextNewIndex = 0;
            int nextRemoveIndex = removeOffset;
            while (nextOriginalIndex < length && nextNewIndex < newLength) {
                while (nextRemoveIndex < removeEnd && remove[nextRemoveIndex++]) { nextOriginalIndex++; } // skip positions to be spliced.
                // Since we make the nextNewIndex < newLength check in the while condition
                // there is no need to include the following break, as is guaranteed not to be true:
                // if (nextOriginalIndex >= length) break; // we reach the final (last positions are to be spliced.
                final int copyStart = nextOriginalIndex;
                while (++nextOriginalIndex < length && (nextRemoveIndex >= removeEnd || !remove[nextRemoveIndex])) { nextRemoveIndex++; }
                final int copyEnd = nextOriginalIndex;
                final int copyLength = copyEnd - copyStart;
                System.arraycopy(source, sourceOffset + copyStart, result, destOffset + nextNewIndex, copyLength);
                nextNewIndex += copyLength;
            }
        }
        return result;
    }

    @SuppressWarnings("unchecked")
    private static <T> T skimArrayBuildResultArray(final T dest, final int destOffset, final Class<T> destClazz, final int newLength) {
        final T result;

        if (dest == null)
            result = (T) Array.newInstance(destClazz.getComponentType(), newLength + destOffset);
        else if (Array.getLength(dest) < newLength + destOffset) {
            result = (T) Array.newInstance(destClazz.getComponentType(),newLength + destOffset);
            if (destOffset > 0) System.arraycopy(dest,0,result,0,destOffset);
        } else
            result = dest;
        return result;
    }

    @SuppressWarnings("unchecked")
    private static <T> Class<T> skimArrayDetermineDestArrayClass(final T dest, Class<T> sourceClazz) {
        final Class<T> destClazz;
        if (dest == null)
            destClazz = sourceClazz;
        else {
            destClazz = (Class<T>) dest.getClass();
            if (destClazz != sourceClazz) {
                if (!destClazz.isArray())
                    throw new IllegalArgumentException("the destination array class must be an array");
                if (sourceClazz.getComponentType().isAssignableFrom(destClazz.getComponentType()))
                    throw new IllegalArgumentException("the provided destination array class cannot contain values from the source due to type incompatibility");
            }
        }
        return destClazz;
    }
}
