package org.broadinstitute.hellbender.utils;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InvalidObjectException;
import java.io.PrintWriter;
import java.io.Writer;
import java.util.Random;

/**
 * Random DNA sequence generator.
 *
 * <p>
 *     Returned bases are always in upper case and one of the valid four nucleotides 'A', 'C', 'G' and 'T'.
 * </p>
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class RandomDNA {

    private final Random random;

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
     * @throws NullPointerException if {@code destination} is {@code null}.
     */
    public void nextBases(final byte[] destination) {
        nextBases(destination, 0, destination.length);
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
        for (int i = offset; i < to; ) {
            int randomInt = random.nextInt();
            for (int j = 0; i < to && j < Integer.SIZE; j += 2) {
                final int ord = randomInt & 0x000000003;
                randomInt >>= 2;
                switch (ord) {
                    case 0:
                        dest[i++] = 'A';
                        break;
                    case 1:
                        dest[i++] = 'C';
                        break;
                    case 2:
                        dest[i++] = 'G';
                        break;
                    case 3:
                        dest[i++] = 'T';
                        break;
                    default:
                        throw new GATKException.ShouldNeverReachHereException("this cannot be happening!!! " + ord);
                }
            }
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

    public File nextFasta(final SAMSequenceDictionary dict, final int basesPerLine)
        throws IOException {
        final File result = File.createTempFile("random-", ".fasta");
        nextFasta(result, dict, basesPerLine);
        return result;
    }

    public void nextFasta(final File out, final SAMSequenceDictionary dict, final int basesPerLine)
            throws IOException {
        final FileWriter writer = new FileWriter(out);
        nextFasta(writer, dict, basesPerLine);
        writer.close();
    }

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
}
