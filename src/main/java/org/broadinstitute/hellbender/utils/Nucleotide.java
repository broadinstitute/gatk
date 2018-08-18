package org.broadinstitute.hellbender.utils;

import javax.validation.constraints.NotNull;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.stream.LongStream;

/**
 * Represents the nucleotide alphabet with support for IUPAC ambiguity codes.
 *
 * <p>
 *    This enumeration not only contains standard (non-ambiguous) nucleotides, but also
 *    values to represent ambiguous and invalid codes.
 * </p>
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public enum Nucleotide {

    // Standard nucleotide codes,
    // and their one-bit-encoding masks CODE(0xTGCA):
    A(0b0001),
    C(0b0010),
    G(0b0100),
    T(0b1000),

    // Extended codes:
    // CODE(included nucs)
    R(A, G), // Purines.
    Y(C, T), // Pyrimidines.
    S(C, G), // Strong nucletoides.
    W(A, T), // Weak nucleotides.
    K(G, T), // Keto nucleotides.
    M(A, C), // Amino nucleotides.
    // The following 4 tri-nucleotide codes don't have a proper long name, they are simply all-except-one.
    B(C, G, T), // Not-A (B follows A)
    D(A, G, T), // Not-C (D follows C)
    H(A, C, T), // Not-G (H follows G)
    V(A, C, G), // Not-V (V follows T)
    // Any
    N(A, C, G, T), // Any/Unknown
    // and X/invalid-call:
    X(); // Invalid.

    // As far as in enum is concern,
    // references to Uracil (U) are considered equivalent to Thymine (T) as they are transcription equivalent.
    public static final Nucleotide U = T;

    // Convenient long form alternative names for some of the enumeration values:

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
    public static final Nucleotide UNKNOWN = N;
    public static final Nucleotide INVALID = X;

    /**
     * List of the standard (non-redundant) nucleotide values in their preferred alphabetical order.
     */
    public static final List<Nucleotide> STANDARD_DNA_BASES = Collections.unmodifiableList(Arrays.asList(A, C, G, T));

    // actually calling values() is costly (creates a new array every time) and often we do just to find out the
    // total number of constants.
    private static final int NUMBER_OF_CONSTANTS;

    private static final Nucleotide[] baseToValue;
    private static final Nucleotide[] maskToValue;

    static {
        final Nucleotide[] values = values();
        NUMBER_OF_CONSTANTS = values.length;
        baseToValue = new Nucleotide[1 << Byte.SIZE];
        maskToValue = new Nucleotide[1 << 4];
        Arrays.fill(baseToValue, INVALID);
        for (final Nucleotide nucleotide : values) {
            baseToValue[nucleotide.lowerCaseByteEncoding & 0xFF]
                    = baseToValue[nucleotide.upperCaseByteEncoding & 0xFF] = nucleotide;
            maskToValue[nucleotide.mask] = nucleotide;
        }
        baseToValue['u'] = baseToValue['U'] = U;
    }

    private final int mask;
    private final boolean isStandard;
    private Nucleotide complement;
    private Nucleotide transition;
    private Nucleotide transversion;

    /**
     * Holds lower-case byte encoding for this nucleotide; {@code 0} for {@link Nucleotide#INVALID}.
     */
    private final byte lowerCaseByteEncoding;

    /**
     * Holds the upper-case byte encoding for this nucleotide; {@code 0} for {@link Nucleotide#INVALID}.
     */
    private final byte upperCaseByteEncoding;

    Nucleotide(final int mask) {
        this.mask = mask;
        isStandard = Integer.bitCount(mask & 0b1111) == 1;
        lowerCaseByteEncoding = (byte) Character.toLowerCase(name().charAt(0));
        upperCaseByteEncoding = (byte) Character.toUpperCase(name().charAt(0));
    }

    Nucleotide(final Nucleotide ... nucs) {
        this(Arrays.stream(nucs).mapToInt(nuc -> nuc.mask).reduce((a, b) -> a | b).orElse(0));
    }

    /**
     * Returns the base that corresponds to this nucleotide.
     * <p>
     *    The base is returned in uppercase.
     * </p>
     * <p>
     *     The {@link #INVALID} nucleotide does not have an actual base then resulting in an exception.
     * </p>
     * @return a valid byte representation for a nucleotide, {@code 0} for {@link Nucleotide#INVALID}.
     */
    public byte encodeAsByte(final boolean upperCase) {
        return upperCase ? upperCaseByteEncoding : lowerCaseByteEncoding;
    }

    /**
     * Returns the nucleotide encoding in a byte using its upper-case representation.
     * @return a valid upper-case byte representation for a nucleotide, {@code 0} for {@link Nucleotide#INVALID}.
     */
    public byte encodeAsByte() {
        return upperCaseByteEncoding;
    }

    /**
     * Returns the nucleotide that corresponds to a particular {@code byte} typed base code.
     * @param base the query base code.
     * @throws IllegalArgumentException if {@code base} is negative.
     * @return never {@code null}, but {@link #INVALID} if the base code does not
     * correspond to a valid nucleotide specification.
     */
    public static Nucleotide decode(final byte base) {
        return baseToValue[Utils.validIndex(base, baseToValue.length)];
    }

    public static Nucleotide decode(final char base) {
        return decode((byte) base);
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

    public boolean isValid() {
        return this != INVALID;
    }

    /**
     * Checks whether this nucleotide code encloses all possible nucleotides for another code.
     * @param other the other nucleotide to compare to.
     * @return {@code true} iff any nucleotide in {@code other} is enclosed in this code.
     */
    public boolean includes(final Nucleotide other) {
        Utils.nonNull(other);
        return other != INVALID && (mask & other.mask) == other.mask;
    }

    public boolean includes(final byte b) {
        return includes(decode(b));
    }

    public Nucleotide intersect(final Nucleotide other) {
        return maskToValue[mask & other.mask];
    }

    /**
     * Checks whether to base encodings make reference to the same {@link #Nucleotide}
     *  instance regardless of their case.
     * <p>
     *     This method is a shorthard for:
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
        return baseToValue[a] == baseToValue[b] && baseToValue[a] != INVALID;
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
        if (complement == null) {
            final int complementMask = ((mask & A.mask) != 0 ? T.mask : 0)
                    | ((mask & T.mask) != 0 ? A.mask : 0)
                    | ((mask & C.mask) != 0 ? G.mask : 0)
                    | ((mask & G.mask) != 0 ? C.mask : 0);
            complement = maskToValue[complementMask];
        }
        return complement;
    }

    /**
     * Returns the complement for a base code.
     * <p>
     *     When an invalid base is provided this method will return the default encoding for the {@link #INVALID} nucleotide.
     * </p>
     * @param b the input base
     * @param upperCase whether to return the uppercase ({@code true}) or the lower case ({@code false}) byte encoding.
     * @return the complement of the input.
     */
    public static byte complement(final byte b, final boolean upperCase) {
        final Nucleotide value = decode(b);
        final Nucleotide compl = value.complement();
        return compl.encodeAsByte(upperCase);
    }

    /**
     * Returns the complement for a base code.
     * <p>
     *     The case of the output will match the case of the input.
     * </p>
     * <p>
     *     When an invalid base is provided this method will return the default encoding for the {@link #INVALID} nucleotide.
     * </p>
     * @param b the input base
     * @return the complement of the input.
     */
    public static byte complement(final byte b) {
        return complement(b, Character.isUpperCase(b));
    }

    /**
     * Returns the instance that would include all possible transition mutations from this one.
     * @return never {@code null}.
     */
    public Nucleotide transition() {
        if (transition == null) {
            final int transitionMask = ((mask & A.mask) != 0 ? G.mask : 0)
                    | ((mask & G.mask) != 0 ? A.mask : 0)
                    | ((mask & C.mask) != 0 ? T.mask : 0)
                    | ((mask & T.mask) != 0 ? C.mask : 0);
            transition = maskToValue[transitionMask];
        }
        return transition;
    }

    /**
     * Returns the instance that would include all possible tranversion mutations from nucleotides included
     * in this one.
     * @return never {@code null}.
     */
    public Nucleotide transversion() {
        if (transversion == null) {
            final int transversionMask = ((mask & PURINE.mask) != 0 ? PYRIMIDINE.mask : 0)
                | ((mask & PYRIMIDINE.mask) != 0 ? PURINE.mask : 0);
            transversion = maskToValue[transversionMask];
        }
        return transversion;
    }

    /**
     * Transvertion mutation toward a strong or a weak base.
     * <p>
     *     This method provides a non-ambiguous alternative to {@link #transversion()} for
     *     concrete nucleotides.
     * </p>
     *
     * @param strong whether the result should be a strong ({@code S: G, C}) or weak ({@code W: A, T}) nucleotide(s).
     * @return nucleotides that may emerged from such a transversion.
     */
    public Nucleotide transversion(final boolean strong) {
        return transversion().intersect(strong ? STRONG : WEAK);
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
            add(decode(base));
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
         * @throws IllegalArgumentException if {@code bases} are null or
         * it contains negative values.
         */
        public final void addAll(final byte ... bases) {
            Utils.nonNull(bases);
            for (final byte base : bases) {
                add(base);
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
            return LongStream.of(counts).sum();
        }
    }
}
