package org.broadinstitute.hellbender.utils;

import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.lang3.StringUtils;
import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.random.Well19937c;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.UserException;

import java.io.File;
import java.lang.reflect.Array;
import java.math.BigInteger;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.*;

public final class Utils {

    private Utils(){}

    /**
     *  Static random number generator and seed.
     */
    private static final long GATK_RANDOM_SEED = 47382911L;
    private static final Random randomGenerator = new Random(GATK_RANDOM_SEED);
    private static final RandomDataGenerator randomDataGenerator = new RandomDataGenerator(new Well19937c(GATK_RANDOM_SEED));

    public static Random getRandomGenerator() { return randomGenerator; }
    public static RandomDataGenerator getRandomDataGenerator() { return randomDataGenerator; }

    public static void resetRandomGenerator() {
        randomGenerator.setSeed(GATK_RANDOM_SEED);
        randomDataGenerator.reSeed(GATK_RANDOM_SEED);
    }

    private static final int TEXT_WARNING_WIDTH = 68;
    private static final String TEXT_WARNING_PREFIX = "* ";
    private static final String TEXT_WARNING_BORDER = StringUtils.repeat('*', TEXT_WARNING_PREFIX.length() + TEXT_WARNING_WIDTH);
    private static final char ESCAPE_CHAR = '\u001B';

    public static final float JAVA_DEFAULT_HASH_LOAD_FACTOR = 0.75f;

    /** our log, which we want to capture anything from this class */
    private static final Logger logger = LogManager.getLogger(Utils.class);

    public static <T> List<T> cons(final T elt, final List<T> l) {
        final List<T> l2 = new ArrayList<>();
        l2.add(elt);
        if (l != null) {
            l2.addAll(l);
        }
        return l2;
    }

    public static void warnUser(final String msg) {
        warnUser(logger, msg);
    }

    public static void warnUser(final Logger logger, final String msg) {
        for (final String line: warnUserLines(msg)) {
            logger.warn(line);
        }
    }

    public static List<String> warnUserLines(final String msg) {
        final List<String> results = new ArrayList<>();
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
                if (space <= 0) {
                    space = TEXT_WARNING_WIDTH;
                }
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
    private static int getLastSpace(final CharSequence message, final int width) {
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
                if (Character.isLetter(c)) {
                    inEscape = false;
                }
            } else if (Character.isWhitespace(c)) {
                lastSpace = currPos;
            }
            currPos++;
        }
        return lastSpace;
    }

    /**
     * Returns a string of the values in an {@link Object} array joined by a separator.
     *
     * @param separator separator character
     * @param objects  the array with values
     *
     * @throws IllegalArgumentException if {@code separator} or {@code objects} is {@code null}.
     * @return a string with the values separated by the separator
     */
    public static String join(final String separator, final Object ... objects) {
        Utils.nonNull(separator, "the separator cannot be null");
        Utils.nonNull(objects, "the value array cannot be null");

        if (objects.length == 0) {
            return "";
        } else {
            final StringBuilder ret = new StringBuilder();
            ret.append(String.valueOf(objects[0]));
            for (int i = 1; i < objects.length; i++) {
                ret.append(separator).append(String.valueOf(objects[i]));
            }
            return ret.toString();
        }
    }

    /**
     * Returns a string of the values in ints joined by separator, such as A,B,C
     *
     * @param separator separator character
     * @param ints   the array with values
     * @return a string with the values separated by the separator
     */
    public static String join(final String separator, final int[] ints) {
        Utils.nonNull(separator, "the separator cannot be null");
        Utils.nonNull(ints, "the ints cannot be null");
        if ( ints.length == 0) {
            return "";
        } else {
            final StringBuilder ret = new StringBuilder();
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
    public static String join(final String separator, final double[] doubles) {
        Utils.nonNull(separator, "the separator cannot be null");
        Utils.nonNull(doubles, "the doubles cannot be null");
        if ( doubles.length == 0) {
            return "";
        } else {
            final StringBuilder ret = new StringBuilder();
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
            {
                return first.toString();
            } else { // full path for 2+ collection that actually need a join
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
        Utils.nonNull(values, "the input array cannot be null");
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
        Utils.nonNull(values, "the input array cannot be null");
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
    @SafeVarargs
    public static <T> List<T> append(final List<T> left, final T ... elts) {
        Utils.nonNull(left, "left is null");
        Utils.nonNull(elts, "the input array cannot be null");
        final List<T> l = new LinkedList<>(left);
        for (final T t : elts){
            Utils.nonNull(t, "t is null");
            l.add(t);
        }
        return l;
    }

    /**
     * Create a new string that's n copies of c
     * @param c the char to duplicate
     * @param nCopies how many copies?
     * @return a string
     */

    public static String dupChar(final char c, final int nCopies) {
        final char[] chars = new char[nCopies];
        Arrays.fill(chars, c);
        return new String(chars);
    }

    /**
     * Create a new byte array that's n copies of b
     * @param b the byte to duplicate
     * @param nCopies how many copies?
     * @return a byte array
     */

    public static byte[] dupBytes(final byte b, final int nCopies) {
        final byte[] bytes = new byte[nCopies];
        Arrays.fill(bytes, b);
        return bytes;
    }


    /**
     * Returns the number of occurrences of a boolean element in a boolean array.
     * @param element
     * @param array cannot be null
     * @return
     */
    public static int countBooleanOccurrences(final boolean element, final boolean[] array) {
        Utils.nonNull(array);
        int count = 0;
        for (final boolean b : array) {
            if (element == b) {
                count++;
            }
        }

        return count;
    }

    /**
     * Splits expressions in command args by spaces and returns the array of expressions.
     * Expressions may use single or double quotes to group any individual expression, but not both.
     * @param args Arguments to parse.
     * @return Parsed expressions.
     */
    public static String[] escapeExpressions(final String args) {
        Utils.nonNull(args);

        // special case for ' and " so we can allow expressions
        if (args.indexOf('\'') != -1) {
            return escapeExpressions(args, "'");
        } else if (args.indexOf('\"') != -1) {
            return escapeExpressions(args, "\"");
        } else {
            return args.trim().split(" +");
        }
    }

    /**
     * Splits expressions in command args by spaces and the supplied delimiter and returns the array of expressions.
     * @param args Arguments to parse.
     * @param delimiter Delimiter for grouping expressions.
     * @return Parsed expressions.
     */
    private static String[] escapeExpressions(final String args, final String delimiter) {
        String[] command = {};
        final String[] split = args.split(delimiter);
        for (int i = 0; i < split.length - 1; i += 2) {
            final String arg = split[i].trim();
            if (arg.length() > 0) { // if the unescaped arg has a size
                command = ArrayUtils.addAll(command, arg.split(" +"));
            }
            command = ArrayUtils.addAll(command, split[i + 1]);
        }
        final String arg = split[split.length - 1].trim();
        if (split.length % 2 == 1 && arg.length() > 0) { // if the last unescaped arg has a size
            command = ArrayUtils.addAll(command, arg.split(" +"));
        }
        return command;
    }


    /**
     * makes an array filled with n copies of the given char.
     */
    public static byte[] repeatChars(final char c, final int n) {
        return repeatBytes((byte)c, n);
    }

    /**
     * makes an array filled with n copies of the given byte.
     */
    public static byte[] repeatBytes(final byte b, final int n) {
        if (n < 0){
            throw new IllegalArgumentException("negative length");
        }
        final byte[] bytes = new byte[n];
        Arrays.fill(bytes, b);
        return bytes;
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
        final List<List<T>> combinations = new ArrayList<>();

        if ( n == 1 ) {
            for ( final T o : objects ) {
                combinations.add(Collections.singletonList(o));
            }
        } else if (n > 1) {
            final List<List<T>> sub = makePermutations(objects, n - 1, withReplacement);
            for ( final List<T> subI : sub ) {
                for ( final T a : objects ) {
                    if ( withReplacement || ! subI.contains(a) ) {
                        combinations.add(Utils.cons(a, subI));
                    }
                }
            }
        }

        return combinations;
    }

    /**
     * @see #calcMD5(byte[])
     */
    public static String calcMD5(final String s) {
        Utils.nonNull(s, "s is null");
        return calcMD5(s.getBytes());
    }

    /**
     * Calculate the md5 for bytes, and return the result as a 32 character string
     *
     * @param bytes the bytes to calculate the md5 of
     * @return the md5 of bytes, as a 32-character long string
     */
    public static String calcMD5(final byte[] bytes) {
        Utils.nonNull(bytes, "the input array cannot be null");
        try {
            final byte[] thedigest = MessageDigest.getInstance("MD5").digest(bytes);
            final BigInteger bigInt = new BigInteger(1, thedigest);

            String md5String = bigInt.toString(16);
            while (md5String.length() < 32) {
                md5String = "0" + md5String; // pad to length 32
            }
            return md5String;
        }
        catch ( final NoSuchAlgorithmException e ) {
            throw new IllegalStateException("MD5 digest algorithm not present");
        }
    }

    /**
     * Checks that an Object {@code object} is not null and returns the same object or throws an {@link IllegalArgumentException}
     * @param object any Object
     * @return the same object
     * @throws IllegalArgumentException if a {@code o == null}
     */
    public static <T> T nonNull(final T object) throws IllegalArgumentException {
        return Utils.nonNull(object, "Null object is not allowed here.");
    }

    /**
     * Checks that an {@link Object} is not {@code null} and returns the same object or throws an {@link IllegalArgumentException}
     * @param object any Object
     * @param message the text message that would be pass to the exception thrown when {@code o == null}.
     * @return the same object
     * @throws IllegalArgumentException if a {@code o == null}
     */
    public static <T> T nonNull(final T object, final String message) throws IllegalArgumentException {
        if (object == null) {
            throw new IllegalArgumentException(message);
        }
        return object;
    }

    /**
     * Checks that the collection does not contain a {@code null} value (throws an {@link IllegalArgumentException} if it does).
     * The implementation calls {@code c.contains(null)} to determine the presence of null.
     * @param c collection
     * @throws IllegalArgumentException if a {@code o == null}
     */
    public static <T> void containsNoNull(final Collection<?> c, final String s) {
        if (c.contains(null)){
            throw new IllegalArgumentException(s);
        }
    }

    /**
     * Checks whether an index is within bounds considering a collection or array of a particular size
     * whose first position index is 0
     * @param index the query index.
     * @param length the collection or array size.
     * @return same value as the input {@code index}.
     */
    public static int validIndex(final int index, final int length) {
        if (index < 0) {
            throw new IllegalArgumentException("the index cannot be negative: " + index);
        } else if (index >= length) {
            throw new IllegalArgumentException("the index points past the last element of the collection or array: " + index + " > " + (length -1));
        }
        return index;
    }

    public static void validateArg(final boolean condition, final String msg){
        if (!condition){
            throw new IllegalArgumentException(msg);
        }
    }

    /**
     * Checks that a user provided file is in fact a regular (i.e. not a directory or a special device) readable file.
     *
     * @param file the input file to test.
     * @throws IllegalArgumentException if {@code file} is {@code null} or {@code argName} is {@code null}.
     * @throws UserException if {@code file} is not a regular file or it cannot be read.
     * @return the same as the input {@code file}.
     */
    public static File regularReadableUserFile(final File file) {
        nonNull(file, "unexpected null file reference");
        if (!file.canRead()) {
            throw new UserException.CouldNotReadInputFile(file.getAbsolutePath(),"the input file does not exist or cannot be read");
        } else if (!file.isFile()) {
            throw new UserException.CouldNotReadInputFile(file.getAbsolutePath(),"the input file is not a regular file");
        } else {
            return file;
        }
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
    public static int optimumHashSize ( final int maxElements ) {
        return (int) (maxElements / JAVA_DEFAULT_HASH_LOAD_FACTOR) + 2;
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
    public static boolean equalRange(final byte[] left, final int leftOffset, final byte[] right, final int rightOffset, final int length) {
        Utils.nonNull(left, "left cannot be null");
        Utils.nonNull(right, "right cannot be null");
        validRange(length, leftOffset, left.length, "left");
        validRange(length, rightOffset, right.length, "right");

        for (int i = 0; i < length; i++) {
            if (left[leftOffset + i] != right[rightOffset + i]) {
                return false;
            }
        }
        return true;
    }

    private static void validRange(final int length, final int offset, final int size, final String msg){
        if (length < 0) {
            throw new IllegalArgumentException(msg + " length cannot be negative");
        }
        if (offset < 0) {
            throw new IllegalArgumentException(msg + " offset cannot be negative");
        }
        if (offset + length > size) {
            throw new IllegalArgumentException(msg + " length goes beyond end of left array");
        }
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
        Utils.nonNull(source, "the source array cannot be null");

        @SuppressWarnings("unchecked")
        final Class<T> sourceClazz = (Class<T>) source.getClass();

        if (!sourceClazz.isArray()) {
            throw new IllegalArgumentException("the source array is not in fact an array instance");
        }
        final int length = Array.getLength(source) - sourceOffset;
        if (length < 0) {
            throw new IllegalArgumentException("the source offset goes beyond the source array length");
        }
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
        Utils.nonNull(source, "the source array cannot be null");
        Utils.nonNull(remove, "the remove array cannot be null");
        if (sourceOffset < 0) {
            throw new IllegalArgumentException("the source array offset cannot be negative");
        }
        if (destOffset < 0) {
            throw new IllegalArgumentException("the destination array offset cannot be negative");
        }
        if (removeOffset < 0) {
            throw new IllegalArgumentException("the remove array offset cannot be negative");
        }
        if (length < 0) {
            throw new IllegalArgumentException("the length provided cannot be negative");
        }

        final int removeLength = Math.min(remove.length - removeOffset,length);

        if (removeLength < 0) {
            throw new IllegalArgumentException("the remove offset provided falls beyond the remove array end");
        }


        @SuppressWarnings("unchecked")
        final Class<T> sourceClazz = (Class<T>) source.getClass();

        if (!sourceClazz.isArray()) {
            throw new IllegalArgumentException("the source array is not in fact an array instance");
        }

        final Class<T> destClazz = skimArrayDetermineDestArrayClass(dest, sourceClazz);

        final int sourceLength = Array.getLength(source);

        if (sourceLength < length + sourceOffset) {
            throw new IllegalArgumentException("the source array is too small considering length and offset");
        }

        // count how many positions are to be removed.

        int removeCount = 0;

        final int removeEnd = removeLength + removeOffset;
        for (int i = removeOffset; i < removeEnd; i++) {
            if (remove[i]) {
                removeCount++;
            }
        }


        final int newLength = length - removeCount;


        @SuppressWarnings("unchecked")
        final T result = skimArrayBuildResultArray(dest, destOffset, destClazz, newLength);
        // No removals, just copy the whole thing.

        if (removeCount == 0) {
            System.arraycopy(source, sourceOffset, result, destOffset, length);
        } else if (length > 0) {  // if length == 0 nothing to do.
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

        if (dest == null) {
            result = (T) Array.newInstance(destClazz.getComponentType(), newLength + destOffset);
        } else if (Array.getLength(dest) < newLength + destOffset) {
            result = (T) Array.newInstance(destClazz.getComponentType(),newLength + destOffset);
            if (destOffset > 0) {
                System.arraycopy(dest, 0, result, 0, destOffset);
            }
        } else {
            result = dest;
        }
        return result;
    }

    @SuppressWarnings("unchecked")
    private static <T> Class<T> skimArrayDetermineDestArrayClass(final T dest, final Class<T> sourceClazz) {
        final Class<T> destClazz;
        if (dest == null) {
            destClazz = sourceClazz;
        } else {
            destClazz = (Class<T>) dest.getClass();
            if (destClazz != sourceClazz) {
                if (!destClazz.isArray()) {
                    throw new IllegalArgumentException("the destination array class must be an array");
                }
                if (sourceClazz.getComponentType().isAssignableFrom(destClazz.getComponentType())) {
                    throw new IllegalArgumentException("the provided destination array class cannot contain values from the source due to type incompatibility");
                }
            }
        }
        return destClazz;
    }

    public static byte [] arrayFromArrayWithLength(final byte[] array, final int length) {
        final byte [] output = new byte[length];
        for (int j = 0; j < length; j++) {
            output[j] = array[(j % array.length)];
        }
        return output;
    }

}
