package org.broadinstitute.hellbender.tools.spark.sv.utils;

import htsjdk.samtools.util.SequenceUtil;
import org.broadinstitute.hellbender.utils.Nucleotide;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;

import java.util.Arrays;

/**
 * Utility class to compose base sequences by appending parts.
 */
public class SequenceBuilder {

    private byte[] bases;

    private int length;

    public SequenceBuilder(final int initialCapacity) {
        length = 0;
        bases = new byte[initialCapacity];
    }

    /**
     * Expands the bases array depending on prospective length expansion.
     * @param lengthExpand number of bases that are going to ab added.
     * @return never {@code null}.
     */
    private SequenceBuilder expandCapacity(final int lengthExpand) {
        if (length + lengthExpand > bases.length) {
            bases = Arrays.copyOf(bases, (length + lengthExpand) << 1);
        }
        return this;
    }

    public SequenceBuilder append(final byte[] source, final int offset, final int length, final boolean reverseComplement) {
        Utils.nonNull(source, "the source array cannot be null");
        ParamUtils.inRange(offset, 0, source.length, "the offset provided is out of range");
        ParamUtils.inRange(length, 0, source.length - offset, "the length provided goes beyond the end of the source array");

        expandCapacity(length);
        System.arraycopy(source, offset, bases, this.length, length);
        if (reverseComplement) {
            SequenceUtil.reverseComplement(bases, this.length, length);
        }
        this.length += length;
        return this;
    }

    public SequenceBuilder append(final byte[] source, final boolean reverseComplement) {
        return append(Utils.nonNull(source), 0, source.length, reverseComplement);
    }

    public SequenceBuilder append(final byte[] source) {
        return append(source, false);
    }

    /**
     * Appends a single base.
     * @param base the base to append.
     * @return this builder.
     */
    public SequenceBuilder append(final byte base) {
        expandCapacity(1);
        bases[length++] = base;
        return this;
    }

    /**
     * Appends bases from a {@link ReferenceBases} instance.
     * @param source the source {@link ReferenceBases}.
     * @param start first position in {@code source} to append (inclusive).
     * @param end last position in {@code source} to append (inclusive).
     * @param reverseComplement whether the appended sequence should be reversed and complemented in the buffer.
     * @return this buffer.
     */
    public SequenceBuilder append(final ReferenceBases source, final int start, final int end, final boolean reverseComplement) {
        Utils.nonNull(source, "the source array cannot be null");
        Utils.validateArg(start <= end, "the start cannot be beyond the end");
        final int expansion = end - start + 1;
        expandCapacity(expansion);
        source.copyBases(bases, length, start, end, false);
        if (reverseComplement) {
            SequenceUtil.reverseComplement(bases, length, expansion);
        }
        length += expansion;
        return this;
    }

    /**
     * Return the current length of the base sequence in the buffer.
     * @return 0 or greater.
     */
    public int getLength() {
        return length;
    }

    /**
     * Creates a copy of the current sequence.
     *
     * @return never {@code null}, an array of exactly {@code getLength()} positions.
     */
    public byte[] toArray() {
        final byte[] result = Arrays.copyOf(bases, length);
        for (int i = 0; i < result.length; i++) {
            if (Nucleotide.decode(result[i]) == Nucleotide.INVALID) {
                throw new IllegalArgumentException("invalid base at " + i + " " + result[i]);
            }
        }
        return result;
    }
}
