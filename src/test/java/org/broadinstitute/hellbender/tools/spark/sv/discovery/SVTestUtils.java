package org.broadinstitute.hellbender.tools.spark.sv.discovery;

import htsjdk.samtools.util.SequenceUtil;

import java.util.Arrays;

public final class SVTestUtils {

    public static byte[] getReverseComplimentCopy(final byte[] sequence) {
        final byte[] sequenceCopy = Arrays.copyOf(sequence, sequence.length);
        SequenceUtil.reverseComplement(sequenceCopy);
        return sequenceCopy;
    }

    public static byte[] makeDummySequence(final int length, byte base) {
        final byte[] result = new byte[length];
        Arrays.fill(result, base);
        return result;
    }
}
