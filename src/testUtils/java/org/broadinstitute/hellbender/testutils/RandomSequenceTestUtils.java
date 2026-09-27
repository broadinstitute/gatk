package org.broadinstitute.hellbender.testutils;

import org.apache.commons.lang3.ArrayUtils;
import org.broadinstitute.hellbender.utils.BaseUtils;

import java.util.Random;

/**
 * Generates random DNA sequences and randomly edited copies of them, for tests that need many related sequence pairs.
 */
public final class RandomSequenceTestUtils {

    private RandomSequenceTestUtils() {}

    /**
     * Returns {@code length} bases drawn uniformly from A, C, G and T.
     */
    public static byte[] randomBases(final Random rng, final int length) {
        final byte[] bases = new byte[length];
        for (int i = 0; i < length; i++) {
            bases[i] = BaseUtils.BASES[rng.nextInt(BaseUtils.BASES.length)];
        }
        return bases;
    }

    /**
     * Returns a copy of {@code bases} with {@code editCount} random edits applied one after another. Six in ten edits
     * are substitutions, two in ten insert between one and {@code maxIndelLength} random bases, and two in ten delete
     * up to {@code maxIndelLength} bases.
     */
    public static byte[] withRandomEdits(final Random rng, final byte[] bases, final int editCount, final int maxIndelLength) {
        byte[] result = bases;
        for (int edit = 0; edit < editCount && result.length > 0; edit++) {
            final int position = rng.nextInt(result.length);
            final int kind = rng.nextInt(10);
            if (kind < 6) {
                result = splice(result, position, position + 1, new byte[]{ substitute(rng, result[position]) });
            } else if (kind < 8) {
                result = splice(result, position, position, randomBases(rng, 1 + rng.nextInt(maxIndelLength)));
            } else {
                final int length = Math.min(1 + rng.nextInt(maxIndelLength), result.length - position);
                result = splice(result, position, position + length, new byte[0]);
            }
        }
        return result;
    }

    /**
     * Returns a base drawn uniformly from the three bases other than {@code base}.
     */
    private static byte substitute(final Random rng, final byte base) {
        final int index = BaseUtils.simpleBaseToBaseIndex(base);
        return BaseUtils.BASES[(index + 1 + rng.nextInt(BaseUtils.BASES.length - 1)) % BaseUtils.BASES.length];
    }

    /**
     * Returns {@code bases} with the range [{@code start}, {@code end}) replaced by {@code replacement}.
     */
    private static byte[] splice(final byte[] bases, final int start, final int end, final byte[] replacement) {
        return ArrayUtils.addAll(ArrayUtils.addAll(ArrayUtils.subarray(bases, 0, start), replacement),
                                 ArrayUtils.subarray(bases, end, bases.length));
    }
}
