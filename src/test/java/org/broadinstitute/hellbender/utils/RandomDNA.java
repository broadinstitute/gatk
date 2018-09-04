package org.broadinstitute.hellbender.utils;

import com.google.api.client.repackaged.com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.util.Random;

/**
 * Random DNA sequence generator.
 *
 * <p>
 *     Returned bases are always in upper case and one of the valid four nucleotides 'A', 'C', 'G' and 'T'.
 * </p>
 *
 * <p>
 *     This class is written in a way that ensure that the sequence of bases returned would only depend on the seed and
 *     not the sequence of calls to different base generating methods.
 * </p>
 * <p>
 *     For example a single call asking for 132 bases using {@link #nextBases(int)} would generate the same sequence
 *     as calling 132 times {@link #nextBase()} if both random instances were initialized with the same seed.
 * </p>
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class RandomDNA {

    /**
     * Number of bases that can be packed in a int (bits in an int divide by 2 bits per base).
     * <p>
     * In practice the code below would work if this constant takes any values between 3 and the actual number
     * of bit pairs in an int (16). The maximum (16) is recommended.
     * </p>
     */
    @VisibleForTesting
    static final int BASES_IN_AN_INT = Integer.SIZE >> 1; // = 16; and is extremely unlikely to ever change.

    /**
     * Convenient constant
     */
    private static final int BASES_IN_AN_INT_MINUS_1 = BASES_IN_AN_INT - 1;

    /**
     * Maximum capacity of the next bytes queue.
     *
     * It must be at least {@link #BASES_IN_AN_INT} and ideally a multiple of that number
     * (in order to avoid unused space in the base-buffer).
     */
    @VisibleForTesting
    static final int NEXT_BASES_MAX_CAPACITY = BASES_IN_AN_INT * 64; // = 1024 bases/bytes  extremelly unlikely to ever change.

    protected final Random random;

    private final byte[] codeToBase = new byte[] {
            Nucleotide.A.encodeAsByte(),
            Nucleotide.C.encodeAsByte(),
            Nucleotide.G.encodeAsByte(),
            Nucleotide.T.encodeAsByte()
    };

    private final byte[] nextBases;
    private int firstBaseIndex;
    private int lastBaseIndex;
    private int nextBasesSize;

    /**
     * Constructs a new random DNA generator.
     * <p>
     * <p>
     * The seed would be the default which would depend on system properties and the current time as
     * described in {@link Random} documentation.
     * </p>
     */
    public RandomDNA() {
        this(new Random());
    }


    /**
     * Creates a new random DNA generator given a random number generator.
     *
     * @param rnd the underlying random number generator.
     * @throws IllegalArgumentException if {@code rnd} is {@code null}.
     */
    public RandomDNA(final Random rnd) {
        random = Utils.nonNull(rnd, "the random number generator cannot be null");
        nextBases = new byte[NEXT_BASES_MAX_CAPACITY];
        firstBaseIndex = 0;
        lastBaseIndex = -1;
        nextBasesSize = 0;
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
     * <p>
     * <p>
     * The whole array will be filled with new base values.
     * </p>
     *
     * @param destination the array to update.
     * @throws IllegalArgumentException if {@code destination} is {@code null}.
     */
    public void nextBases(final byte[] destination) {
        nextBases(Utils.nonNull(destination, "the destination array must not be null"), 0, destination.length);
    }

    /**
     * Puts a random sequence of base into an byte array.
     *
     * @param dest   the destination array.
     * @param offset first position on the destination array to fill.
     * @param length how many based to put into the array.
     * @throws IllegalArgumentException if {@code dest} is {@code null} or the pair {@code offset} and {@code length}
     *                                  do not specify a valid index sub-interval in {@code dest}.
     */
    public void nextBases(final byte[] dest, final int offset, final int length) {
        Utils.nonNull(dest, "input array cannot be null");
        ParamUtils.inRange(offset, 0, dest.length, "the offset provided is not a valid index given the destination size");
        ParamUtils.isPositiveOrZero(length, "the input length must be positive");
        final int to = offset + length;
        ParamUtils.inRange(to, 0, dest.length, "the length provided goes beyond the end of the destination array");
        int remaining = length;
        int nextDestIndex = offset;
        while (remaining > 0) {
            final int available = loadNextBases(remaining);
            for (int i = 0; i < available; i++) {
                dest[nextDestIndex++] = removeNextBase();
            }
            remaining -= available;
        }
    }

    /**
     * Load new bases to the next bases queue to try to satisfy a required number of next-bases.
     *
     * <p>
     *     It might be that there is not enough space in the buffer to satisfy the target, in which case this will return
     * something close to the capacity of the buffer, expecting the caller to clear the buffer and repeat the request
     * for the remainder.
     * </p>
     *
     * @param targetSize the number of bases that the caller would like to have at its disposal.
     * @return the amount of bases available given the current content of the next-bases queue, its max capacity
     * and the requested target.
     */
    private int loadNextBases(final int targetSize) {
        final int currentNumberOfBases = nextBasesSize;
        if (currentNumberOfBases >= targetSize) {
            return targetSize;
        } else if (NEXT_BASES_MAX_CAPACITY - currentNumberOfBases < BASES_IN_AN_INT) {
            // If we try to load more bases it would go over the max capacity.
            return currentNumberOfBases;
        } else  {
            final int result = Math.min(NEXT_BASES_MAX_CAPACITY, targetSize);
            final int increaseInInts = (result - currentNumberOfBases + BASES_IN_AN_INT_MINUS_1) / BASES_IN_AN_INT;
            for (int i = 0; i < increaseInInts; i++) {
                int randomInt = random.nextInt();
                addNextBase(codeToBase[randomInt & 0x03]);
                int j = 1;
                do {
                    randomInt >>>= 2;
                    addNextBase(codeToBase[randomInt & 0x03]);
                } while(++j < BASES_IN_AN_INT_MINUS_1);
                randomInt >>>= 2;
                addNextBase(codeToBase[randomInt]); // for the last base there is no need to bother with the mask.
            }
            return result;
        }
    }

    /**
     * Returns a random DNA sequence of bases.
     *
     * @param size the length of the sequence.
     * @return never {@code null}.
     * @throws IllegalArgumentException if {@code size} is negative.
     */
    public byte[] nextBases(final int size) {
        Utils.validateArg(size >= 0, "the size cannot be negative");
        final byte[] result = new byte[size];
        nextBases(result);
        return result;
    }

    /**
     * Create a random reference and writes in FASTA format into a temporal file.
     * <p>
     *     The output file is instantiated using {@link File#createTempFile}.
     * </p>
     * <p>
     *     The invoking code is responsible to manage and dispose of the output file.
     * </p>
     * @param dict the dictionary indicating the number of contigs and their lengths.
     * @param basesPerLine number of base to print in each line of the output FASTA file.
     * @return the temporal file with the random reference.
     * @throws IOException if such an exception was thrown while accessing and writing into the temporal file.
     * @throws IllegalArgumentException if {@code dict} is {@code null} or {@code basesPerLine} is 0 or negative.
     */
    public File nextFasta(final SAMSequenceDictionary dict, final int basesPerLine)
        throws IOException {
        final File result = File.createTempFile("random-", ".fasta");
        nextFasta(result, dict, basesPerLine);
        return result;
    }

    /**
     * Creates a random reference and writes it in FASTA format into a file.
     * @param out the output file.
     * @param dict the dictionary indicating the number of contigs and their lengths.
     * @param basesPerLine number of base to print in each line of the output FASTA file.
     *
     * @throws IOException if such an exception was thrown while accessing and writing into the temporal file.
     * @throws IllegalArgumentException if {@code dict} is {@code null}, or {@code out } is {@code null}
     *    or {@code basesPerLine} is 0 or negative.
     */
    public void nextFasta(final File out, final SAMSequenceDictionary dict, final int basesPerLine)
            throws IOException {
        Utils.nonNull(out);
        try (final FileWriter writer = new FileWriter(out)) {
            nextFasta(writer, dict, basesPerLine);
        }
    }

    /**
     * Creates a random reference and writes it in FASTA format into a {@link Writer}.
     * @param out the output writer.
     * @param dict the dictionary indicating the number of contigs and their lengths.
     * @param basesPerLine number of base to print in each line of the output FASTA file.
     *
     * @throws IOException if such an exception was thrown while accessing and writing into the temporal file.
     * @throws IllegalArgumentException if {@code dict} is {@code null}, or {@code out } is {@code null}
     *    or {@code basesPerLine} is 0 or negative.
     */
    public void nextFasta(final Writer out, final SAMSequenceDictionary dict, final int basesPerLine)
            throws IOException {
        Utils.nonNull(out);
        Utils.nonNull(dict);
        ParamUtils.isPositive(basesPerLine, "number of base per line must be strictly positive: " + basesPerLine);
        final byte[] buffer = new byte[basesPerLine];
        final String lineSeparator = System.lineSeparator();
        for (final SAMSequenceRecord sequence : dict.getSequences()) {
            int pendingBases = sequence.getSequenceLength();
            out.append(">").append(sequence.getSequenceName()).append(lineSeparator);
            while (pendingBases > 0) {
                final int lineLength = pendingBases < basesPerLine ? pendingBases : basesPerLine;
                nextBases(buffer, 0, lineLength);
                out.append(new String(buffer, 0, lineLength)).append(lineSeparator);
                pendingBases -= lineLength;
            }
        }
    }

    public byte nextBase() {
        loadNextBases(1); // this never will fail to have 1 available base.
        return removeNextBase();
    }


    private void addNextBase(final byte base) {
        lastBaseIndex++;
        if (lastBaseIndex == NEXT_BASES_MAX_CAPACITY) {
            lastBaseIndex = 0;
        }
        nextBases[lastBaseIndex] = base;
        nextBasesSize++;
    }

    private byte removeNextBase() {
        final byte result = nextBases[firstBaseIndex++];
        if (firstBaseIndex == NEXT_BASES_MAX_CAPACITY) {
            firstBaseIndex = 0;
        }
        nextBasesSize--;
        return result;
    }
}
