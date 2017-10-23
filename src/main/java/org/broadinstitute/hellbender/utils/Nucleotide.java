package org.broadinstitute.hellbender.utils;

import java.util.Arrays;
import java.util.stream.LongStream;

/**
 * Represents the nucleotide alphabet.
 *
 * <p>
 *    This enumeration not only contains concrete nucleotides, but also
 *    values to represent ambiguous and invalid codes.
 * </p>
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public enum Nucleotide {
    A(true, false, false, false),
    C(false, true, false, false),
    G(false, false, true, false),
    T(false, false, false, true),
    R(true, false, true, false),
    Y(false, true, false, true),
    S(false, true, true, false),
    W(true, false, false, true),
    K(false, false, true, true),
    M(true, false, true, false),
    B(false, true, true, true),
    D(true, false, true, true),
    H(true, true, false, true),
    V(true, true, true, false),
    N(true, true, true, true),
    INVALID(false, false, false, false);

    public static final Nucleotide U = T;
    public static final Nucleotide X = N;

    private static final Nucleotide[] baseToValue = new Nucleotide[Byte.MAX_VALUE + 1];

    static {
        Arrays.fill(baseToValue, INVALID);
        for (final Nucleotide nucleotide : values()) {
            baseToValue[nucleotide.lowerCaseByteEncoding] = baseToValue[nucleotide.upperCaseByteEncoding] = nucleotide;
        }
        baseToValue['u'] = baseToValue['U'] = U;
        baseToValue['x'] = baseToValue['X'] = X;
    }

    private static final int A_MASK = 1;
    private static final int C_MASK = A_MASK << 1;
    private static final int G_MASK = C_MASK << 1;
    private static final int T_MASK = G_MASK << 1;

    private final int acgtMask;
    private final boolean isConcrete;

    /**
     * Holds lower-case byte encoding for this nucleotide; {@code 0} for {@link Nucleotide#INVALID}.
     */
    private final byte lowerCaseByteEncoding;

    /**
     * Holds the upper-case byte encoding for this nucleotide; {@code 0} for {@link Nucleotide#INVALID}.
     */
    private final byte upperCaseByteEncoding;

    Nucleotide(final boolean a, final boolean c, final boolean g, final boolean t) {
        acgtMask = (a ? A_MASK : 0) | (c ? C_MASK : 0) | (g ? G_MASK : 0) | (t ? T_MASK : 0);
        isConcrete = acgtMask == A_MASK || acgtMask == C_MASK || acgtMask == G_MASK || acgtMask == T_MASK;
        lowerCaseByteEncoding = acgtMask == 0 ? (byte) 0 : (byte) name().charAt(0);
        upperCaseByteEncoding = acgtMask == 0 ? (byte) 0 : (byte) Character.toUpperCase(name().charAt(0));
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
    public byte encodeAsByte(final boolean lowerCase) {
        return lowerCase ? lowerCaseByteEncoding : upperCaseByteEncoding;
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

    /**
     * Checks whether the nucleotide refer to a concrete (rather than ambiguous) base.
     * @return {@code true} iff this is a concrete nucleotide.
     */
    public boolean isConcrete() {
        return isConcrete;
    }

    /**
     * Checks whether the nucleotide refer to an ambiguous base.
     * @return {@code true} iff this is an ambiguous nucleotide.
     */
    public boolean isAmbiguous() {
        return !isConcrete && this != INVALID;
    }

    public boolean includes(final Nucleotide other) {
        Utils.nonNull(other);
        if (this == INVALID || other == INVALID) {
            return false;
        } else {
            return ((this.acgtMask & other.acgtMask) == other.acgtMask);
        }
    }

    /**
     * Helper class to count the number of occurrences of each nucleotide in
     * a sequence.
     */
    public static class Counter {
        private final long[] counts;

        /**
         * Creates a new counter with all counts set to 0.
         */
        public Counter() {
            counts = new long[Nucleotide.values().length];
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
        public void addAll(final byte[] bases) {
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

        public long sum() {
            return LongStream.of(counts).sum();
        }
    }
}
