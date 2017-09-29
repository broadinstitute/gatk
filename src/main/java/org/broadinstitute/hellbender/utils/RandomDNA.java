package org.broadinstitute.hellbender.utils;

import java.util.Random;

/**
 * Random DNA sequence generator.
 *
 * <p>
 *     Returned bases are always in upper case and one of the valid four nocleotides 'A', 'C', 'G' and 'T'.
 * </p>
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class RandomDNA {

    private final Random random;

    /**
     * Constructs a new random DNA generator.
     *
     * <p>
     *     The seed would be the default which would depend on system properties and the current time as
     *     described in {@link Random} documentation.
     * </p>
     */
    public RandomDNA() {
        this(new Random());
    }


    /**
     * Creates a new random DNA generator given a random number generator.
     * @param rnd the underlying random number generator.
     *
     * @throws IllegalArgumentException if {@code rnd} is {@code null}.
     */
    public RandomDNA(final Random rnd) {
        random = Utils.nonNull(rnd,"the random number generator cannot be null");
    }

    /**
     * Constructs a new random DNA generator providing a seed.
     *
     * @param seed the random number generator seed.
     */
    public RandomDNA(final long seed) {
        this(new Random(seed));
    }

    /**
     * Updates the content of a byte array with a random base sequence.
     *
     * <p>
     *     The whole array will be filled with new base values.
     * </p>
     *
     * @param destination the array to update.
     *
     * @throws NullPointerException if {@code destination} is {@code null}.
     */
    public void nextBases(final byte[] destination) {
        random.nextBytes(destination);
        for (int i = 0; i < destination.length; i++) {
            final int ord = destination[i] & 0x03;
            switch (ord) {
                case 0: destination[i] = 'A'; break;
                case 1: destination[i] = 'C'; break;
                case 2: destination[i] = 'G'; break;
                case 3: destination[i] = 'T'; break;
                default: throw new IllegalStateException("this cannot be happening!!!");
            }
        }
    }

    /**
     * Returns a random RNA sequence of bases.
     * @param size the length of the sequence.
     *
     * @throws IllegalArgumentException if {@code size} is negative.
     * @return never {@code null}.
     */
    public byte[] nextBases(final int size) {
        Utils.validateArg(size >= 0, "the size cannot be negative");
        final byte[] result = new byte[size];
        nextBases(result);
        return result;
    }


}
