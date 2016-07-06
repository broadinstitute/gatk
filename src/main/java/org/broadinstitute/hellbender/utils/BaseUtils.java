package org.broadinstitute.hellbender.utils;

import org.broadinstitute.hellbender.exceptions.UserException;

import java.util.Arrays;
import java.util.Comparator;
import java.util.Random;

/**
 * BaseUtils contains some basic utilities for manipulating nucleotides.
 */
public final class BaseUtils {

    public static boolean isNBase(byte base) {
        return base == 'N' || base == 'n';
    }

    public enum Base {
        A ('A'),
        C ('C'),
        G ('G'),
        T ('T'),
        N ('N'),
        D ('D');

        public byte base;

        private Base(final char base) {
            this.base = (byte)base;
        }
    }

    // todo -- add this to the generalized base abstraction using the Base enum.
    public static final byte[] BASES = {'A', 'C', 'G', 'T'};
    public static final char[] BASE_CHARS = {'A', 'C', 'G', 'T'};

    private static final int[] baseIndexMap = new int[256];
    static {
        Arrays.fill(baseIndexMap, -1);
        baseIndexMap['A'] = Base.A.ordinal();
        baseIndexMap['a'] = Base.A.ordinal();
        baseIndexMap['*'] = Base.A.ordinal();    // the wildcard character counts as an A
        baseIndexMap['C'] = Base.C.ordinal();
        baseIndexMap['c'] = Base.C.ordinal();
        baseIndexMap['G'] = Base.G.ordinal();
        baseIndexMap['g'] = Base.G.ordinal();
        baseIndexMap['T'] = Base.T.ordinal();
        baseIndexMap['t'] = Base.T.ordinal();
    }

    private static final int[] baseIndexWithIupacMap = Arrays.copyOf(baseIndexMap, baseIndexMap.length);
    static {
        baseIndexWithIupacMap['*'] = -1;    // the wildcard character is bad
        baseIndexWithIupacMap['N'] = Base.N.ordinal();
        baseIndexWithIupacMap['n'] = Base.N.ordinal();
        baseIndexWithIupacMap['R'] = Base.N.ordinal();
        baseIndexWithIupacMap['r'] = Base.N.ordinal();
        baseIndexWithIupacMap['Y'] = Base.N.ordinal();
        baseIndexWithIupacMap['y'] = Base.N.ordinal();
        baseIndexWithIupacMap['M'] = Base.N.ordinal();
        baseIndexWithIupacMap['m'] = Base.N.ordinal();
        baseIndexWithIupacMap['K'] = Base.N.ordinal();
        baseIndexWithIupacMap['k'] = Base.N.ordinal();
        baseIndexWithIupacMap['W'] = Base.N.ordinal();
        baseIndexWithIupacMap['w'] = Base.N.ordinal();
        baseIndexWithIupacMap['S'] = Base.N.ordinal();
        baseIndexWithIupacMap['s'] = Base.N.ordinal();
        baseIndexWithIupacMap['B'] = Base.N.ordinal();
        baseIndexWithIupacMap['b'] = Base.N.ordinal();
        baseIndexWithIupacMap['D'] = Base.N.ordinal();
        baseIndexWithIupacMap['d'] = Base.N.ordinal();
        baseIndexWithIupacMap['H'] = Base.N.ordinal();
        baseIndexWithIupacMap['h'] = Base.N.ordinal();
        baseIndexWithIupacMap['V'] = Base.N.ordinal();
        baseIndexWithIupacMap['v'] = Base.N.ordinal();
    }

    /// In genetics, a transition is a mutation changing a purine to another purine nucleotide (A <-> G) or
    // a pyrimidine to another pyrimidine nucleotide (C <-> T).
    // Approximately two out of every three single nucleotide polymorphisms (SNPs) are transitions.
    public enum BaseSubstitutionType {
        TRANSITION,         // A <-> G or C <-> T
        TRANSVERSION
    }

    /**
     * Returns the base substitution type of the 2 state SNP
     *
     * @param base1
     * @param base2
     * @return
     */
    public static BaseSubstitutionType SNPSubstitutionType(final byte base1, final byte base2) {
        final BaseSubstitutionType t = isTransition(base1, base2) ? BaseSubstitutionType.TRANSITION : BaseSubstitutionType.TRANSVERSION;
        //System.out.printf("SNPSubstitutionType( char %c, char %c ) => %s%n", base1, base2, t);
        return t;
    }

    public static boolean isTransition(final byte base1, final byte base2) {
        final int b1 = simpleBaseToBaseIndex(base1);
        final int b2 = simpleBaseToBaseIndex(base2);
        return b1 == Base.A.ordinal() && b2 == Base.G.ordinal() || b1 == Base.G.ordinal() && b2 == Base.A.ordinal() ||
                b1 == Base.C.ordinal() && b2 == Base.T.ordinal() || b1 == Base.T.ordinal() && b2 == Base.C.ordinal();
    }

    /**
     * Private constructor.  No instantiating this class!
     */
    private BaseUtils() {}

    public static boolean basesAreEqual(final byte base1, final byte base2) {
        return simpleBaseToBaseIndex(base1) == simpleBaseToBaseIndex(base2);
    }

    public static byte[] convertIUPACtoN(final byte[] bases, final boolean errorOnBadReferenceBase, final boolean ignoreConversionOfFirstByte) {
        final int length = bases.length;
        final int start = ignoreConversionOfFirstByte ? 1 : 0;

        for ( int i = start; i < length; i++ ) {
            final int baseIndex = baseIndexWithIupacMap[bases[i]];
            if ( baseIndex == Base.N.ordinal() ) {
                bases[i] = 'N';
            } else if ( errorOnBadReferenceBase && baseIndex == -1 ) {
                throw new UserException.BadInput("We encountered a non-standard non-IUPAC base in the provided reference: '" + bases[i] + "'");
            }
        }
        return bases;
    }

    /**
     * Converts a simple base to a base index
     *
     * @param base [AaCcGgTt]
     * @return 0, 1, 2, 3, or -1 if the base can't be understood
     */
    public static int simpleBaseToBaseIndex(final byte base) {
        Utils.validateArg( base >= 0 && base < 256,
                "Non-standard bases were encountered in either the input reference or BAM file(s)");
        return baseIndexMap[base];
    }

    /**
     * Returns true iff the base represented by the byte is a 'regular' base (ACGT or *).
     */
    public static boolean isRegularBase( final byte base ) {
        return simpleBaseToBaseIndex(base) != -1;
    }

    /**
     * Returns true iff all bases are 'regular' {@link #isRegularBase}.
     */
    public static boolean isAllRegularBases( final byte[] bases ) {
        for( final byte base : bases) {
            if( !isRegularBase(base) ) { return false; }
        }
        return true;
    }

    /**
     * Converts a base index to a simple base
     *
     * @param baseIndex 0, 1, 2, 3
     * @return A, C, G, T, or '.' if the index can't be understood
     */
    public static byte baseIndexToSimpleBase(final int baseIndex) {
        switch (baseIndex) {
            case 0:
                return 'A';
            case 1:
                return 'C';
            case 2:
                return 'G';
            case 3:
                return 'T';
            default:
                return '.';
        }
    }

    /**
     * Return the complement (A <-> T or C <-> G) of a base, or the specified base if it can't be complemented (i.e. an ambiguous base).
     *
     * @param base the base [AaCcGgTt]
     * @return the complementary base, or the input base if it's not one of the understood ones
     */
    public static byte simpleComplement(final byte base) {
        switch (base) {
            case 'A':
            case 'a':
                return 'T';
            case 'C':
            case 'c':
                return 'G';
            case 'G':
            case 'g':
                return 'C';
            case 'T':
            case 't':
                return 'A';
            default:
                return base;
        }
    }

    /**
     * Reverse complement a byte array of bases (that is, chars casted to bytes, *not* base indices in byte form)
     *
     * @param bases the byte array of bases
     * @return the reverse complement of the base byte array
     */
    public static byte[] simpleReverseComplement(final byte[] bases) {
        final byte[] rcbases = new byte[bases.length];

        for (int i = 0; i < bases.length; i++) {
            rcbases[i] = simpleComplement(bases[bases.length - 1 - i]);
        }

        return rcbases;
    }

    // --------------------------------------------------------------------------------
    //
    // random bases
    //
    // --------------------------------------------------------------------------------


    /**
     * Fill an array section with random bases.
     *
     * @param dest array to fill.
     * @param fromIndex first index to be filled (inclusive).
     * @param toIndex index after last to be filled (exclusive).
     *
     * @throws IllegalArgumentException if {@code dest} is {@code null},
     *              {@code fromIndex} or {@code toIndex} is negative,
     *              {@code fromIndex} or {@code toIndex} are greater than {@code dest} length,
     *              or {@code fromIndex} greater than {@code toIndex}.
     */
    public static void fillWithRandomBases(final byte[] dest, final int fromIndex, final int toIndex) {
        final Random rnd = Utils.getRandomGenerator();
        Utils.nonNull(dest, "the dest array cannot be null");
        Utils.validateArg(fromIndex <= toIndex, "fromIndex cannot be larger than toIndex");
        Utils.validateArg(fromIndex >= 0, "both indexes must be positive");
        Utils.validateArg(toIndex <= dest.length, "both indexes must be less or equal to the destination array length");

        new IndexRange(fromIndex, toIndex).forEach(i -> dest[i] = baseIndexToSimpleBase(rnd.nextInt(4)));
    }

    public static byte getComplement(final byte base) {
        switch(base) {
            case 'a':
            case 'A':
                return 'T';
            case 'c':
            case 'C':
                return 'G';
            case 'g':
            case 'G':
                return 'C';
            case 't':
            case 'T':
                return 'A';
            case 'n':
            case 'N':
                return 'N';
            default:
                throw new IllegalArgumentException("base must be A, C, G or T. " + (char) base + " is not a valid base.");
        }
    }

    /**
     * Lexicographical sorting of base arrays {@link Comparator}.
     */
    public static final Comparator<byte[]> BASES_COMPARATOR = (byte[] o1, byte[] o2) ->
    {
            Utils.nonNull(o1, "o1");
            Utils.nonNull(o2, "o2");
            final int minLength = Math.min(o1.length,o2.length);
            for (int i = 0; i < minLength; i++) {
                final int cmp = Byte.compare(o1[i],o2[i]);
                if (cmp != 0) {
                    return cmp;
                }
            }
            return Integer.compare(o1.length, o2.length);
    };
}
