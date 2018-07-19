package org.broadinstitute.hellbender.utils;

import com.google.common.collect.AbstractIterator;
import com.google.common.collect.Iterators;
import com.google.common.primitives.Ints;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.tribble.util.ParsingUtils;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.lang3.StringUtils;
import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.random.Well19937c;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.exceptions.GATKException;

import javax.annotation.Nullable;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.lang.reflect.Array;
import java.math.BigInteger;
import java.nio.file.Files;
import java.nio.file.Path;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.time.ZonedDateTime;
import java.time.format.DateTimeFormatter;
import java.time.format.FormatStyle;
import java.util.*;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.function.Function;
import java.util.function.Supplier;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;

public final class Utils {

    /**
     * Comparator for strings that sorts null first;
     */
    public static final Comparator<? super String> COMPARE_STRINGS_NULLS_FIRST = Comparator.nullsFirst(Comparator.naturalOrder());

    private Utils(){}

    private final static DateTimeFormatter longDateTimeFormatter = DateTimeFormatter.ofLocalizedDateTime(FormatStyle.LONG);

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
     * Create a new string thats a n duplicate copies of s
     * @param s the string to duplicate
     * @param nCopies how many copies?
     * @return a string
     */
    public static String dupString(final String s, int nCopies) {
        if ( s == null || s.equals("") ) { throw new IllegalArgumentException("Bad s " + s); }
        if ( nCopies < 0 ) { throw new IllegalArgumentException("nCopies must be >= 0 but got " + nCopies); }

        final StringBuilder b = new StringBuilder();
        for ( int i = 0; i < nCopies; i++ ) {
            b.append(s);
        }
        return b.toString();
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
            if (!arg.isEmpty()) { // if the unescaped arg has a size
                command = ArrayUtils.addAll(command, arg.split(" +"));
            }
            command = ArrayUtils.addAll(command, split[i + 1]);
        }
        final String arg = split[split.length - 1].trim();
        if (split.length % 2 == 1 && !arg.isEmpty()) { // if the last unescaped arg has a size
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
            throw new IllegalStateException("MD5 digest algorithm not present", e);
        }
    }

    /**
     * Calculates the MD5 for the specified file and returns it as a String
     *
     * Warning: this loads the whole file into memory, so it's not suitable
     * for large files.
     *
     * @param file file whose MD5 to calculate
     * @return file's MD5 in String form
     * @throws IOException if the file could not be read
     */
    public static String calculateFileMD5( final File file ) throws IOException{
        return calculatePathMD5(file.toPath());
    }

    /**
     * Calculates the MD5 for the specified file and returns it as a String
     *
     * Warning: this loads the whole file into memory, so it's not suitable
     * for large files.
     *
     * @param path file whose MD5 to calculate
     * @return file's MD5 in String form
     * @throws IOException if the file could not be read
     */
    public static String calculatePathMD5(final Path path) throws IOException{
        // This doesn't have as nice error messages as FileUtils, but it's close.
        String fname = path.toUri().toString();
        if (!Files.exists(path)) {
            throw new FileNotFoundException("File '" + fname + "' does not exist");
        }
        if (Files.isDirectory(path)) {
            throw new IOException("File '" + fname + "' exists but is a directory");
        }
        if (!Files.isRegularFile(path)) {
            throw new IOException("File '" + fname + "' exists but is not a regular file");
        }
        return Utils.calcMD5(Files.readAllBytes(path));
    }

    /**
     * Checks that an Object {@code object} is not null and returns the same object or throws an {@link IllegalArgumentException}
     * @param object any Object
     * @return the same object
     * @throws IllegalArgumentException if a {@code o == null}
     */
    public static <T> T nonNull(final T object) {
        return Utils.nonNull(object, "Null object is not allowed here.");
    }

    /**
     * Checks that an {@link Object} is not {@code null} and returns the same object or throws an {@link IllegalArgumentException}
     * @param object any Object
     * @param message the text message that would be passed to the exception thrown when {@code o == null}.
     * @return the same object
     * @throws IllegalArgumentException if a {@code o == null}
     */
    public static <T> T nonNull(final T object, final String message) {
        if (object == null) {
            throw new IllegalArgumentException(message);
        }
        return object;
    }

    /**
     * Checks that an {@link Object} is not {@code null} and returns the same object or throws an {@link IllegalArgumentException}
     * @param object any Object
     * @param message the text message that would be passed to the exception thrown when {@code o == null}.
     * @return the same object
     * @throws IllegalArgumentException if a {@code o == null}
     */
    public static <T> T nonNull(final T object, final Supplier<String> message) {
        if (object == null) {
            throw new IllegalArgumentException(message.get());
        }
        return object;
    }


    /**
     * Checks that a {@link Collection} is not {@code null} and that it is not empty.
     * If it's non-null and non-empty it returns the input, otherwise it throws an {@link IllegalArgumentException}
     * @param collection any Collection
     * @param message a message to include in the output
     * @return the original collection
     * @throws IllegalArgumentException if collection is null or empty
     */
    public static <I, T extends Collection<I>> T nonEmpty(T collection, String message){
        nonNull(collection, "The collection is null: " + message);
        if(collection.isEmpty()){
            throw new IllegalArgumentException("The collection is empty: " + message);
        } else {
            return collection;
        }
    }

    /**
     * Checks that a {@link Collection} is not {@code null} and that it is not empty.
     * If it's non-null and non-empty it returns the true
     * @param collection any Collection
     * @return true if the collection exists and has elements
     */
    public static boolean isNonEmpty(Collection<?> collection){
        return collection != null && !collection.isEmpty();
    }

    /**
     * Checks that a {@link String} is not {@code null} and that it is not empty.
     * If it's non-null and non-empty it returns the input, otherwise it throws an {@link IllegalArgumentException}
     * @param string any String
     * @param message a message to include in the output
     * @return the original string
     * @throws IllegalArgumentException if string is null or empty
     */
    public static String nonEmpty(String string, String message){
        nonNull(string, "The string is null: " + message);
        if(string.isEmpty()){
            throw new IllegalArgumentException("The string is empty: " + message);
        } else {
            return string;
        }
    }

    /**
     * Checks that a {@link String} is not {@code null} and that it is not empty.
     * If it's non-null and non-empty it returns the input, otherwise it throws an {@link IllegalArgumentException}
     * @param string any String
     * @return the original string
     * @throws IllegalArgumentException if string is null or empty
     */
    public static String nonEmpty(final String string){
        return nonEmpty(string, "string must not be null or empty");
    }

    /**
     * Checks that a {@link Collection} is not {@code null} and that it is not empty.
     * If it's non-null and non-empty it returns the input, otherwise it throws an {@link IllegalArgumentException}
     * @param collection any Collection
     * @return the original collection
     * @throws IllegalArgumentException if collection is null or empty
     */
    public static  <I, T extends Collection<I>> T nonEmpty(T collection){
        return nonEmpty(collection, "collection must not be null or empty.");
    }

    /**
     * Checks that the collection does not contain a {@code null} value (throws an {@link IllegalArgumentException} if it does).
     * @param collection collection
     * @param message the text message that would be pass to the exception thrown when c contains a null.
     * @throws IllegalArgumentException if collection is null or contains any null elements
     */
    public static void containsNoNull(final Collection<?> collection, final String message) {
        Utils.nonNull(collection, message);
        //cannot use Collection.contains(null) here because this throws a NullPointerException when used with many Sets
        if (collection.stream().anyMatch(v -> v == null)){
            throw new IllegalArgumentException(message);
        }
    }

    /**
     * Checks that the collection does not contain a duplicate value (throws an {@link IllegalArgumentException} if it does).
     * The implementation creates a {@link Set} as an intermediate step or detecting duplicates and returns this Set because
     * it is sometimes useful to do so.
     *
     * @param c collection
     * @param message A message to emit in case of error, in addition to reporting the first duplicate value found.
     * @throws IllegalArgumentException if a {@code o == null}
     */
    public static <E> Set<E> checkForDuplicatesAndReturnSet(final Collection<E> c, final String message) {
        final Set<E> set = new LinkedHashSet<>();
        for (final E element : c) {
            if (!set.add(element)) {
                throw new IllegalArgumentException(String.format(message + "  Value %s appears more than once.", element.toString()));
            }
        }
        return set;
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

    public static void validateArg(final boolean condition, final Supplier<String> msg){
        if (!condition){
            throw new IllegalArgumentException(msg.get());
        }
    }

    /**
     * Check a condition that should always be true and throw an {@link IllegalStateException} if false.  If msg is not a
     * String literal i.e. if it requires computation, use the Supplier<String> version, below.
     */
    public static void validate(final boolean condition, final String msg){
        if (!condition){
            throw new IllegalStateException(msg);
        }
    }

    /**
     * Check a condition that should always be true and throw an {@link IllegalStateException} if false.
     */
    public static void validate(final boolean condition, final Supplier<String> msg){
        if (!condition){
            throw new IllegalStateException(msg.get());
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
     * Skims out positions of an array returning a shorter one with the remaining positions in the same order.
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
     * @param destOffset the first position where to copy the skimmed array values.
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

    /**
     * Checks if the read header contains any reads groups from non-Illumina and issue a warning of that's the case.
     */
    public static void warnOnNonIlluminaReadGroups(final SAMFileHeader readsHeader, final Logger logger) {
        Utils.nonNull(readsHeader, "header");
        Utils.nonNull(logger, "logger");
        if (readsHeader.getReadGroups().stream().anyMatch(rg -> NGSPlatform.fromReadGroupPL(rg.getPlatform()) !=  NGSPlatform.ILLUMINA)){
            logger.warn("This tool has only been well tested on ILLUMINA-based sequencing data. For other data use at your own risk.");
        }
    }

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
     * Find the last occurrence of the query sequence in the reference sequence
     *
     * Returns the index of the last occurrence or -1 if the query sequence is not found
     *
     * @param reference the reference sequence
     * @param query the query sequence
     */
    public static int lastIndexOf(final byte[] reference, final byte[] query) {
        int queryLength = query.length;

        // start search from the last possible matching position and search to the left
        for (int r = reference.length - queryLength; r >= 0; r--) {
            int q = 0;
            while (q < queryLength && reference[r+q] == query[q]) {
                q++;
            }
            if (q == queryLength) {
                return r;
            }
        }
        return -1;
    }

    /**
     * Simple wrapper for sticking elements of a int[] array into a List<Integer>
     * @param ar - the array whose elements should be listified
     * @return - a List<Integer> where each element has the same value as the corresponding index in @ar
     */
    public static List<Integer> listFromPrimitives(final int[] ar) {
        Utils.nonNull(ar);
        return Ints.asList(ar);
    }

    /**
     * Concatenates a series of {@link Iterator}s (all of the same type) into a single {@link Iterator}.
     * @param iterator an {@link Iterator} of {@link Iterator}s
     * @param <T> the type of the iterator
     * @return an {@link Iterator} over the underlying {@link Iterator}s
     */
    public static <T> Iterator<T> concatIterators(final Iterator<? extends Iterable<T>> iterator) {
        Utils.nonNull(iterator, "iterator");

        return new AbstractIterator<T>() {
            Iterator<T> subIterator;
            @Override
            protected T computeNext() {
                if (subIterator != null && subIterator.hasNext()) {
                    return subIterator.next();
                }
                while (iterator.hasNext()) {
                    subIterator = iterator.next().iterator();
                    if (subIterator.hasNext()) {
                        return subIterator.next();
                    }
                }
                return endOfData();
            }
        };
    }

    public static <T> Stream<T> stream(final Iterable<T> iterable) {
        return StreamSupport.stream(iterable.spliterator(), false);
    }

    public static <T> Stream<T> stream(final Iterator<T> iterator) {
        return stream(() -> iterator);
    }

    /**
     * Like Guava's {@link Iterators#transform(Iterator, com.google.common.base.Function)}, but runs a fixed number
     * ({@code numThreads}) of transformations in parallel, while maintaining ordering of the output iterator.
     * This is useful if the transformations are CPU intensive.
     */
    public static <F, T> Iterator<T> transformParallel(final Iterator<F> fromIterator, final Function<F, T> function, final int numThreads) {
        Utils.nonNull(fromIterator, "fromIterator");
        Utils.nonNull(function, "function");
        Utils.validateArg(numThreads >= 1, "numThreads must be at least 1");

        if (numThreads == 1) { // defer to Guava for single-threaded case
            return Iterators.transform(fromIterator, new com.google.common.base.Function<F, T>() {
                @Nullable
                @Override
                public T apply(@Nullable final F input) {
                    return function.apply(input);
                }
            });
        }
        // use an executor service for the multi-threaded case
        final ExecutorService executorService = Executors.newFixedThreadPool(numThreads);
        final Queue<Future<T>> futures = new LinkedList<>();
        return new AbstractIterator<T>() {
            @Override
            protected T computeNext() {
                try {
                    while (fromIterator.hasNext()) {
                        if (futures.size() == numThreads) {
                            return futures.remove().get();
                        }
                        final F next = fromIterator.next();
                        final Future<T> future = executorService.submit(() -> function.apply(next));
                        futures.add(future);
                    }
                    if (!futures.isEmpty()) {
                        return futures.remove().get();
                    }
                    executorService.shutdown();
                    return endOfData();
                } catch (InterruptedException | ExecutionException e) {
                    throw new GATKException("Problem running task", e);
                }
            }
        };
    }

    /** Gets duplicated items in the collection. */
    public static <T> Set<T> getDuplicatedItems(final Collection<T> objects) {
        final Set<T> unique = new HashSet<>();
        return objects.stream()
                .filter(name -> !unique.add(name))
                .collect(Collectors.toSet());
    }

    /**
     * Return the given {@code dateTime} formatted as string for display.
     * @param dateTime the date/time to be formatted
     * @return String representing the {@code dateTime}.
     */
    public static String getDateTimeForDisplay(final ZonedDateTime dateTime) {
        return dateTime.format(longDateTimeFormatter);
    }

    /**
     * Set the Locale to US English so that numbers will always be formatted in the US style.
     */
    public static void forceJVMLocaleToUSEnglish() {
        Locale.setDefault(Locale.US);
    }

    /**
     * Streams and sorts a collection of objects and returns the integer median entry of the sorted list
     * @param values List of sortable entries from which to select the median
     */
    public static <T extends Comparable<?>> T getMedianValue(List<T> values) {
        final List<T> sorted = values.stream().sorted().collect(Collectors.toList());
        return sorted.get(sorted.size() / 2);
    }


    /**
     * Splits a String using indexOf instead of regex to speed things up.
     * This method produces the same results as {@link String#split(String)} and {@code String.split(String, 0)},
     * but has been measured to be ~2x faster (see {@code StringSplitSpeedUnitTest} for details).
     *
     * @param str       the string to split.
     * @param delimiter the delimiter used to split the string.
     * @return A {@link List} of {@link String} tokens.
     */
    public static List<String> split(final String str, final char delimiter) {

        final List<String> tokens;

        if ( str.isEmpty() ) {
            tokens = new ArrayList<>(1);
            tokens.add("");
        }
        else {
            tokens = ParsingUtils.split(str, delimiter);
            removeTrailingEmptyStringsFromEnd(tokens);
        }

        return tokens;
    }

    /**
     * Splits a String using indexOf instead of regex to speed things up.
     * If given an empty delimiter, will return each character in the string as a token.
     * This method produces the same results as {@link String#split(String)} and {@code String.split(String, 0)},
     * but has been measured to be ~2x faster (see {@code StringSplitSpeedUnitTest} for details).
     *
     * @param str       the string to split.
     * @param delimiter the delimiter used to split the string.
     * @return A {@link List} of {@link String} tokens.
     */
    public static List<String> split(final String str, final String delimiter) {
        // This is 10 because the ArrayList default capacity is 10 (but private).
        return split(str, delimiter, 10);
    }

    /**
     * Splits a given {@link String} using {@link String#indexOf} instead of regex to speed things up.
     * If given an empty delimiter, will return each character in the string as a token.
     * This method produces the same results as {@link String#split(String)} and {@code String.split(String, 0)},
     * but has been measured to be ~2x faster (see {@code StringSplitSpeedUnitTest} for details).
     *
     * @param str               The {@link String} to split.
     * @param delimiter         The delimiter used to split the {@link String}.
     * @param expectedNumTokens The number of tokens expected (used to initialize the capacity of the {@link ArrayList}).
     * @return A {@link List} of {@link String} tokens.
     */
    private static List<String> split(final String str, final String delimiter, final int expectedNumTokens) {
        final List<String> result;

        if ( str.isEmpty() ) {
            result = new ArrayList<>(1);
            result.add("");
        }
        else if ( delimiter.isEmpty() ) {
            result = new ArrayList<>(str.length());
            for ( int i = 0; i < str.length(); ++i ) {
                result.add(str.substring(i, i + 1));
            }
        }
        else if ( delimiter.length() == 1 ) {
            result = split(str, delimiter.charAt(0));
        }
        else {
            result = new ArrayList<>(expectedNumTokens);

            int delimiterIdx = -1;
            int tokenStartIdx = delimiterIdx + 1;
            do {
                delimiterIdx = str.indexOf(delimiter, tokenStartIdx);
                final String token = (delimiterIdx != -1 ? str.substring(tokenStartIdx, delimiterIdx) : str.substring(tokenStartIdx));
                result.add(token);
                tokenStartIdx = delimiterIdx + delimiter.length();
            } while ( delimiterIdx != -1 );

            removeTrailingEmptyStringsFromEnd(result);
        }

        return result;
    }

    private static void removeTrailingEmptyStringsFromEnd(final List<String> result) {
        // Remove all trailing empty strings to emulate the behavior of String.split:
        // We remove items from the end of the list to our index
        // so that we can take advantage of better performance of removing items from the end
        // of certain concrete lists:
        while ( (!result.isEmpty()) && (result.get(result.size() - 1).isEmpty()) ) {
            result.remove(result.size() - 1);
        }
    }

    /**
     * Take a map of a value to a list and reverse it.  Note that no assumptions of uniqueness are made, so returned
     *  values are also lists.
     *
     *  <p>For example:</p>
     *
     *  Input:<br/>
     *  k ->  {a,b} <br/>
     *  j ->  {a} <br/>
     *
     *  Output:<br/>
     *  a ->  {k,j} <br/>
     *  b ->  {k} <br/>
     *
     *  Any sorting in the input map will be lost in the output.
     *
     * @param somethingToListMap a map from a value to a list of values.  Never {@code null}
     * @param <T> class of the key of the input
     * @param <U> class of the values in the list of the input
     * @return A new mapping from class of values to set of keys.  Never {@code null}
     */
    public static <T,U> Map<U, Set<T>> getReverseValueToListMap(final Map<T, List<U>> somethingToListMap) {
        final Map<U, Set<T>> result = new HashMap<>();

        for (final Map.Entry<T, List<U>> entry : somethingToListMap.entrySet()) {
            entry.getValue().forEach(v -> result.computeIfAbsent(v, k -> new HashSet<>()).add(entry.getKey()));
        }

        return result;
    }

    /**
     * Convenience function that formats a percentage as a %.2f string
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
     * Given a collection of strings and a collection of regular expressions, generates the set of strings that match
     * any expression
     * @param sourceValues collection of strings from which to to select
     * @param filterExpressions list of expressions to use for matching
     * @param exactMatch If true match filters exactly, otherwise use as both exact and regular expressions
     * @return A new set strings from sourceValues that satisfy at least one of the expressions in sampleExpressions
     */
    public static Set<String> filterCollectionByExpressions(final Collection<String> sourceValues, final Collection<String> filterExpressions, final boolean exactMatch) {
        Utils.nonNull(filterExpressions);
        Utils.nonNull(sourceValues);

        final Set<String> filteredValues = new LinkedHashSet<>();

        Collection<Pattern> patterns = null;
        if (!exactMatch) {
            patterns = compilePatterns(filterExpressions);
        }
        for (final String value : sourceValues) {
            if (filterExpressions.contains(value)) {
                filteredValues.add(value);
            } else if (!exactMatch) {
                for (final Pattern pattern : patterns) {
                    if (pattern.matcher(value).find()) {
                        filteredValues.add(value);
                        break;
                    }
                }
            }
        }

        return filteredValues;
    }

    private static Collection<Pattern> compilePatterns(final Collection<String> filters) {
        final Collection<Pattern> patterns = new ArrayList<Pattern>();
        for (final String filter: filters) {
            patterns.add(Pattern.compile(filter));
        }
        return patterns;
    }
}
