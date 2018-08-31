package org.broadinstitute.hellbender.utils;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * Represents the nucleotide alphabet with support for IUPAC ambiguity codes.
 *
 * <p>
 *    This enumeration not only contains standard (non-ambiguous) nucleotides, but also
 *    contains ambiguous nucleotides, as well as a code {@link #X} (a.k.a. {@link #INVALID})
 *    for invalid nucleotide calls.
 * </p>
 *
 * <p>
 *     You can query whether a value refers to a non-ambiguous nucleotide with {@link #isStandard()} or
 *     {@link #isAmbiguous()} whichever is most convenient. Notice that the special value {@link #X}
 *     is neither of those.
 * </p>
 *
 * <p>
 *     Querying the {@link #X} value for its {@link #complement}, {@link #transition} or
 *     {@link #transversion} or using it in other operations
 *     such as {@link #intersect} will return {@link #X}; similar to {@link Double#NaN} in
 *     {@code double} arithmetic.
 * </p>
 *
 * <p>
 *     For naming consistency it is recommended to use {@link #decode} and {@link #encodeAsByte}
 *     or {@link #encodeAsString} methods to translate byte/char and string encodings from and
 *     into values of this enum over the inherited {@link #toString}, {@link #name} or {@link #valueOf}.
 * </p>
 *
 * <p>
 *     Although the canonical names for values use the single letter IUPAC
 *     encodings, this class provides convenient longer form names constant aliases
 *     (e.g. {@link #ADENINE} for {@link #A}, {@link #PURINE} for {@link #R}, etc.).
 * </p>
 * <p>
 *     Uracil and Thymine are considered equivalent in this enum with {@link #T} as the canonical name.
 * </p>
 * <p>
 *     Finally, notice that there is no code of the "gap nucleotide" that may appear in aligned sequences as in fact
 *     that is not a nucleotide. A base encoding using the typical gap representation such as '.' or '-' would
 *     be interpreted as an {@link #INVALID} (i.e. {@link #X}) call which is probably not what you want.
 *     So code to support those will need to do so outside this {@code enum}.
 * </p>
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public enum Nucleotide {

    // Standard nucleotide codes,
    // and their one-bit-encoding masks CODE(0bTGCA):
    A(0b0001),
    C(0b0010),
    G(0b0100),
    T(0b1000),

    // Extended codes:
    // CODE(included nucs)
    R(A, G), // Purines.
    Y(C, T), // Pyrimidines.
    S(C, G), // Strong nucleotides.
    W(A, T), // Weak nucleotides.
    K(G, T), // Keto nucleotides.
    M(A, C), // Amino nucleotides.
    // The following 4 tri-nucleotide codes don't have a proper long name, they are simply "all-except-one"
    // codes:
    B(C, G, T), // Not-A (B follows A)
    D(A, G, T), // Not-C (D follows C)
    H(A, C, T), // Not-G (H follows G)
    V(A, C, G), // Not-V (V follows T)
    // Any
    N(A, C, G, T), // Any/Unknown

    // And X/invalid-call:
    X();

    // As far as this enum is concern,
    // references to Uracil (U) are considered equivalent to Thymine (T) as they are transcription equivalent.
    // nucleotides.
    public static final Nucleotide U = T;

    // Convenient constants with long form alternative names for some of the enumeration values:

    // Long form standard nucleotide names.
    public static final Nucleotide ADENINE = A;
    public static final Nucleotide CYTOSINE = C;
    public static final Nucleotide GUANINE = G;
    public static final Nucleotide THYMINE = T;
    public static final Nucleotide URACIL = U;

    // Ambiguous nucleotide groups with proper long form names:
    public static final Nucleotide STRONG = S;
    public static final Nucleotide WEAK = W;
    public static final Nucleotide PURINE = R;
    public static final Nucleotide PYRIMIDINE = Y;
    public static final Nucleotide AMINO = M;
    public static final Nucleotide KETO = K;
    public static final Nucleotide ANY = N;
    public static final Nucleotide INVALID = X;

    /**
     * List of the standard (non-redundant) nucleotide values in their preferred alphabetical order.
     */
    public static final List<Nucleotide> STANDARD_BASES = Collections.unmodifiableList(Arrays.asList(A, C, G, T));

    // Since calling values() is costly (creates a new array every time) and often we do it just to find out the
    // total number of constants is best to cache it in a constant.
    private static final int NUMBER_OF_CONSTANTS;

    /**
     * Values indexed by their unsigned byte encodings. Non-valid encodings point to {@link #INVALID}.
     */
    private static final Nucleotide[] baseToValue;

    /**
     * Values indexed by their mask.
     */
    private static final Nucleotide[] maskToValue;

    /**
     * Value ordinal indexed by their unsigned byte ecodings. Non-valid encodings point to {@link #INVALID}
     * (thru its ordinal).
     */
    private static final int[] baseToOrdinal;

    static {
        final Nucleotide[] values = values();
        NUMBER_OF_CONSTANTS = values.length;
        baseToValue = new Nucleotide[1 << Byte.SIZE];
        maskToValue = new Nucleotide[1 << STANDARD_BASES.size()];
        baseToOrdinal = new int[1 << Byte.SIZE];
        Arrays.fill(baseToValue, INVALID);
        Arrays.fill(baseToOrdinal, INVALID.ordinal());
        for (final Nucleotide nucleotide : values) {
            // Notice that {@code "x & 0xFF"} is needed instead of {@code "(int)x" as
            // we want the unsigned value (e.g. 255 rather than -1).
            // This is repeated through this class code.
            final int lowerCaseIndex = nucleotide.lowerCaseByteEncoding & 0xFF;
            final int upperCaseIndex = nucleotide.upperCaseByteEncoding & 0xFF;
            maskToValue[nucleotide.mask] = nucleotide;
            baseToValue[lowerCaseIndex] = baseToValue[upperCaseIndex] = nucleotide;
            baseToOrdinal[lowerCaseIndex] = baseToOrdinal[upperCaseIndex] = nucleotide.ordinal();
        }
        // need to do u and U here as they are just aliases to T.
        baseToValue['u' & 0xFF] = baseToValue['U' & 0xFF] = U;
        baseToOrdinal['u' & 0xFF] = baseToOrdinal['U' & 0xFF] = U.ordinal();
    }

    private final int mask;
    private final boolean isStandard;

    // Some properties initialized after construction as these depend on some static arrays
    // defined above.
    private Nucleotide complement;
    private Nucleotide transition;
    private Nucleotide transversion;

    static {
        for (final Nucleotide value : values()) {
            value.finalizeInitialization();
        }
    }

    /**
     * Holds the lower-case byte encoding for this nucleotide.
     * This is typically the lower-case version of the enum constant name.
     */
    private final byte lowerCaseByteEncoding;

    /**
     * Holds the lower-case {@link String} representation for this nucleotide.
     * This is typically the lower-case version of the enum constant name.
     */
    private final String lowerCaseStringEncoding;

    /**
     * Holds the lower-case {@code char} representation for this nucleotide.
     * This is typically the lower-case version of the enum constant name only character.
     */
    private final char lowerCaseCharEncoding;

    /**
     * Holds the upper-case byte encoding for this nucleotide.
     * This is typically the upper-case version of the enum constant name.
     */
    private final byte upperCaseByteEncoding;

    /**
     * Holds the upper-case {@code char} representation for this nucleotide.
     * This is typically the upper-case version of the enum constant name only character.
     */
    private final char upperCaseCharEncoding;

    /**
     * Construct a nucleotide given its mask.
     * @param mask the mask.
     */
    Nucleotide(final int mask) {
        this.mask = mask;
        isStandard = Integer.bitCount(mask & 0b1111) == 1;
        lowerCaseByteEncoding = (byte) Character.toLowerCase(name().charAt(0));
        lowerCaseCharEncoding = Character.toLowerCase(name().charAt(0));
        upperCaseByteEncoding = (byte) Character.toUpperCase(name().charAt(0));
        upperCaseCharEncoding = Character.toUpperCase(name().charAt(0));
        lowerCaseStringEncoding = name().toLowerCase();
    }

    /**
     * Construct a nucleotide given the other codes that it would include.
     * @param nucs the nucleotides to include.
     */
    Nucleotide(final Nucleotide ... nucs) {
        this(Arrays.stream(nucs).mapToInt(nuc -> nuc.mask).reduce((a, b) -> a | b).orElse(0));
    }

    /**
     * Returns the {@code byte} typed encoding that corresponds to this nucleotide.
     * @param upperCase whether to return the upper- or lower-case {@code byte} representation.
     * @return a valid and exclusive {@code byte} representation for a nucleotide.
     */
    public byte encodeAsByte(final boolean upperCase) {
        return upperCase ? upperCaseByteEncoding : lowerCaseByteEncoding;
    }

    /**
     * Returns the {@code char} typed encoding that corresponds to this nucleotide.
     * @param upperCase whether to return the upper- or lower-case {@code char} representation.
     * @return a valid and exclusive {@code char} representation for a nucleotide.
     */
    public char encodeAsChar(final boolean upperCase) {
        return upperCase ? upperCaseCharEncoding : lowerCaseCharEncoding;
    }

    /**
     * Returns this nucleotide's exclusive upper-case {@code byte} encoding.
     * @return <i>ditto</i>.
     */
    public byte encodeAsByte() {
        return upperCaseByteEncoding;
    }

    /**
     * Returns the nucleotide's exclusive upper-case {@code char} encoding.
     * @return <i>ditto</i>.
     */
    public char encodeAsChar() {
        return upperCaseCharEncoding;
    }

    /**
     * Returns the nucleotide's exclusive upper-case {@code String} encoding.
     * @return <i>ditto</i>.
     */
    public String encodeAsString() {
        return toString();
    }

    /**
     * Returns the nucleotide's exclusive {@link String} typed encoding.
     * @param upperCase whether the upper or lower-case representation should be returned.
     * @return a valid and exclusive {@link String} representation for this nucleotide.
     */
    public String encodeAsString(final boolean upperCase) {
        return upperCase ? toString() : lowerCaseStringEncoding;
    }

    /**
     * Returns the nucleotide that corresponds to a particular {@code byte} typed base code.
     * @param base the query base code.
     * @return never {@code null}, but {@link #INVALID} if the base code does not
     * correspond to a valid nucleotide specification.
     */
    public static Nucleotide decode(final byte base) {
        return baseToValue[base & 0xFF];
    }

    /**
     * Returns the nucleotide that corresponds to a particular {@code char} typed base code.
     * @param ch the query base code.
     * @return never {@code null}, but {@link #INVALID} if the base code does not correspond
     * to a valid nucleotide specification.
     */
    public static Nucleotide decode(final char ch) {
        if ((ch & 0xFF00) != 0) {
            return INVALID;
        } else {
            return baseToValue[ch & 0xFF];
        }
    }

    /**
     * Transform a single-letter character string into the corresponding nucleotide.
     * <p>
     *    {@code Null}, empty or multi-letter input will result in an {@link IllegalArgumentException}.
     *    These are not simply invalid encodings as the fact that are not a single character is
     *    an indication of a probable bug.
     * </p>
     *
     * @param seq the input character sequence to transform into.
     * @return never {@code null}, perhaps {@link #INVALID} to indicate that the input is not a valid
     * single letter encoding encoding.
     */
    public static Nucleotide decode(final CharSequence seq) {
        Utils.nonNull(seq, "the input character sequence must not be null");
        if (seq.length() != 1) {
            throw new IllegalArgumentException("the input character sequence must be exactly one character long");
        } else {
            return decode(seq.charAt(0));
        }
    }

    /**
     * Checks whether the nucleotide refers to a concrete (rather than ambiguous) base.
     * @return {@code true} iff this is a concrete nucleotide.
     */
    public boolean isStandard() {
        return isStandard;
    }

    /**
     * Checks whether the nucleotide refer to an ambiguous base.
     * @return {@code true} iff this is an ambiguous nucleotide.
     */
    public boolean isAmbiguous() {
        return !isStandard && this != INVALID;
    }

    /**
     * Whether this nucleotide code is valid or not.
     * @return {@code true} iff valid.
     */
    public boolean isValid() {
        return this != INVALID;
    }

    /**
     * Checks whether this nucleotide code encloses all possible nucleotides for another code.
     * @param other the other nucleotide to compare to.
     * @return {@code true} iff every nucleotide in {@code other} is enclosed in this code.
     */
    public boolean includes(final Nucleotide other) {
        Utils.nonNull(other);
        return other != INVALID && (mask & other.mask) == other.mask;
    }

    /**
     * Checks whether this nucleotide code encloses all possible nucleotides for another code.
     * @param b the other nucleotide to compare to encoded as a byte.
     * @return {@code true} iff every nucleotide in {@code other} is enclosed in this code.
     */
    public boolean includes(final byte b) {
        return includes(decode(b));
    }

    /**
     * Returns the nucleotide code that include all and only the nucleotides that are
     * included by this another code.
     * @param other the other nucleotide code.
     * @throws IllegalArgumentException if {@code other} is {@code null}.
     * @return never {@code null}. Returns {@link #INVALID} if the intersection does not contain
     * any nucleotide.
     */
    public Nucleotide intersect(final Nucleotide other) {
        Utils.nonNull(other, "the other nucleotide cannot be null");
        return maskToValue[mask & other.mask];
    }

    /**
     * Checks whether two nucleotides intersect given their byte encodings.
     * @param a first nucleotide.
     * @param b second nucleotide.
     * @return {@code true} iff the input nucleotides intersect.
     */
    public static boolean intersect(final byte a, final byte b) {
        return (baseToValue[0xFF & a].mask & baseToValue[0xFF & b].mask) != 0;
    }

    /**
     * Checks whether two base encodings make reference to the same {@link #Nucleotide}
     *  instance regardless of their case.
     * <p>
     *     This method is a shorthand for:
     *     <pre>{@link #decode}(a){@link #same(Nucleotide) same}({@link #decode}(b)) </pre>.
     * </p>
     *
     *  <p>
     *      The order of the inputs is not relevant, therefore {@code same(a, b) == same(b, a)} for any
     *      given {@code a} and {@code b}.
     *  </p>
     *  <p>
     *      Notice that if either or both input bases make reference to an invalid nucleotide (i.e. <pre> {@link #decode}(x) == {@link #INVALID}},
     *      this method will return {@code false} even if {@code a == b}.
     *  </p>
     * @param a the first base to compare (however order is not relevant).
     * @param b the second base to compare (however order is not relevant).
     * @return {@code true} iff {@code {@link #decode}}.same({@link #decode}(b))}}
     */
    public static boolean same(final byte a, final byte b) {
        return baseToValue[a & 0xFF] == baseToValue[b & 0xFF] && baseToValue[a & 0xFF] != INVALID;
    }

    /**
     * Checks whether this and another {@link #Nucleotide} make reference to the same nucleotide(s).
     * <p>
     *     In contrast with {@link #equals}, this method will return {@code false} if any of the two, this
     *     or the input nucleotide is the {@link #INVALID} enum value. So even <pre>{@link #INVALID}.same({@link #INVALID})</pre>
     *     will return {@code false}.
     * </p>
     *
     * @param other the other nucleotide.
     * @return {@code true} iff this and the input nucleotide make reference to the same nucleotides.
     */
    public boolean same(final Nucleotide other) {
        return this == other && this != INVALID;
    }

    /**
     * Returns the complement nucleotide code for this one.
     * <p>
     *     For ambiguous nucleotide codes, this will return the ambiguous code that encloses the complement of
     *     each possible nucleotide in this code.
     * </p>
     * <p>
     *     The complement of the {@link #INVALID} nucleotide is itself.
     * </p>
     * @return never {@code null}.
     */
    public Nucleotide complement() {
        return complement;
    }

    /**
     * Returns the complement for a base code.
     * <p>
     *     When an invalid base is provided this method will return the input byte (lower- or upper-cased depending on that
     *     flag value).
     * </p>
     * @param b the input base
     * @param upperCase whether to return the uppercase ({@code true}) or the lower case ({@code false}) byte encoding.
     * @return the complement of the input.
     */
    public static byte complement(final byte b, final boolean upperCase) {
        final Nucleotide compl = baseToValue[b & 0xFF].complement;
        return compl != INVALID
                ? (upperCase ? compl.upperCaseByteEncoding : compl.lowerCaseByteEncoding)
                : (byte) ( upperCase ? Character.toUpperCase(b) : Character.toLowerCase(b));
    }

    /**
     * Returns the complement for a base code.
     * <p>
     *     The case of the output will match the case of the input.
     * </p>
     * <p>
     *     When an invalid base is provided this method will return the input base byte.
     * </p>
     * @param b the input base
     * @return the complement of the input.
     */
    public static byte complement(final byte b) {
        final Nucleotide compl = baseToValue[b & 0xFF].complement;
        return compl != INVALID
                ? (Character.isUpperCase(b) ? compl.upperCaseByteEncoding : compl.lowerCaseByteEncoding)
                : b;
    }

    /**
     * Returns the instance that would include all possible transition mutations from this one.
     * @return never {@code null}.
     */
    public Nucleotide transition() {
        return transition;
    }

    /**
     * Returns the instance that would include all possible tranversion mutations from nucleotides included
     * in this one.
     * @return never {@code null}.
     */
    public Nucleotide transversion() {
        return transversion;
    }

    /**
     * Calculate and set the complement, transition and transversion using the #maskToValue array.
     */
    private void finalizeInitialization() {
       // set the complement.
       final int complementMask = ((mask & A.mask) != 0 ? T.mask : 0)
                        | ((mask & T.mask) != 0 ? A.mask : 0)
                        | ((mask & C.mask) != 0 ? G.mask : 0)
                        | ((mask & G.mask) != 0 ? C.mask : 0);
       complement = maskToValue[complementMask];
       // set the transversion.
       final int transversionMask = ((mask & PURINE.mask) != 0 ? PYRIMIDINE.mask : 0)
                    | ((mask & PYRIMIDINE.mask) != 0 ? PURINE.mask : 0);
       transversion = maskToValue[transversionMask];
       // set the transition.
       final int transitionMask = ((mask & A.mask) != 0 ? G.mask : 0)
                    | ((mask & G.mask) != 0 ? A.mask : 0)
                    | ((mask & C.mask) != 0 ? T.mask : 0)
                    | ((mask & T.mask) != 0 ? C.mask : 0);
       transition = maskToValue[transitionMask];
    }

    /**
     * Transversion mutation toward a strong or a weak base.
     * <p>
     *     This method provides a non-ambiguous alternative to {@link #transversion()} for
     *     concrete nucleotides.
     * </p>
     *
     * @param strong whether the result should be a strong ({@code S: G, C}) or weak ({@code W: A, T}) nucleotide(s).
     * @return nucleotides that may emerged from such a transversion.
     */
    public Nucleotide transversion(final boolean strong) {
        return transversion.intersect(strong ? STRONG : WEAK);
    }

    /**
     * Helper class to count the number of occurrences of each nucleotide code in
     * a sequence.
     */
    public static final class Counter {

        private final long[] counts;

        /**
         * Creates a new counter with all counts set to 0.
         */
        public Counter() {
            counts = new long[NUMBER_OF_CONSTANTS];
        }

        /**
         * Increases by 1 the count for a nucleotide.
         * @param nucleotide the target nucleotide.
         * @throws IllegalArgumentException if nucleotide is {@code null}.
         */
        public void add(final Nucleotide nucleotide) {
            counts[Utils.nonNull(nucleotide).ordinal()]++;
        }

        /**
         * Increases the nucleotide that corresponds to the input base own count by 1.
         * @param base the base code.
         * @throws IllegalArgumentException if {@code base} is {@code negative}.
         */
        public void add(final byte base) {
            counts[baseToOrdinal[base & 0xFF]]++;
        }

        public void add(final char base) {
            if ((base & 0xFF00) != 0) {
                counts[INVALID.ordinal()]++;
            } else {
                counts[baseToOrdinal[base & 0xFF]]++;
            }
        }

        /**
         * Returns the current count for a given nucleotide.
         * @param nucleotide the query nucleotide.
         * @throws IllegalArgumentException if {@code nucleotide} is {@code null}.
         * @return 0 or greater.
         */
        public long get(final Nucleotide nucleotide) {
            return counts[Utils.nonNull(nucleotide).ordinal()];
        }

        /**
         * Increase by one the count for a nucleotide for each
         * occurrence of such in the input byte array base codes.
         * @param bases the input base codes.
         * @throws IllegalArgumentException if {@code bases} is null.
         */
        public final void addAll(final byte ... bases) {
            Utils.nonNull(bases);
            for (final byte base : bases) {
                counts[baseToOrdinal[base & 0xFF]]++;
            }
        }

        /**
         * Increase by one the count for a nucleotide for each
         * occurrence of such in the input char array base codes.
         * @param bases the input base codes.
         * @throws IllegalArgumentException if {@code bases} is null.
         */
        public final void addAll(final char ... bases) {
            Utils.nonNull(bases);
            for (final char base : bases) {
                if ((base & 0xFF00) != 0) {
                    counts[INVALID.ordinal()]++;
                } else {
                    counts[baseToOrdinal[base & 0xFF]]++;
                }
            }
        }

        /**
         * Increase by one the count for a nucleotide for each
         * occurrence of such in the input {@link CharSequence}.
         * @param bases the input bases sequence.
         * @throws IllegalArgumentException if the input is {@code null}.
         */
        public final void addAll(final CharSequence bases) {
            Utils.nonNull(bases);
            for (int i = 0; i < bases.length(); i++) {
                final char base = bases.charAt(i);
                if ((base & 0xFF00) != 0) {
                    counts[INVALID.ordinal()]++;
                } else {
                    counts[baseToOrdinal[base & 0xFF]]++;
                }
            }
        }

        /**
         * Reset all the counts to 0.
         */
        public void clear() {
            Arrays.fill(counts, 0);
        }

        /**
         * Return the total count of all nucleotide constants.
         * @return 0 or greater.
         */
        public long sum() {
            return MathUtils.sum(counts);
        }
    }
}
