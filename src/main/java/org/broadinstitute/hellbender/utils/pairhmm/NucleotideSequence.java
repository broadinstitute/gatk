package org.broadinstitute.hellbender.utils.pairhmm;

import org.broadinstitute.hellbender.utils.Nucleotide;

public interface NucleotideSequence {

    long length();

    default Nucleotide nucleotideAt(final long position) {
        return Nucleotide.decode(byteAt(position));
    }

    byte byteAt(final long position);

    default String subString(final long bestBeg, final long bestEnd) {
        final byte[] bytes = new byte[(int) (bestEnd - bestBeg + 1)];
        for (int i = 0, offset = 0; i < bytes.length; i++, offset++) {
            bytes[i] = byteAt(bestBeg + offset);
        }
        return new String(bytes);
    }

    default int copyBytesAt(byte[] dest, int destOffset, long location, final int length) {
        final long seqLength = this.length();
        int i, j;
        long k;
        for (i = 0, j = destOffset, k = location; i < length && k <= seqLength; i++, k++, j++) {
            dest[j] = byteAt(k);
        }
        return i;
    }

    default byte[] bytesAt(long location, final int length) {
        final byte[] result = new byte[length];
        if (copyBytesAt(result, 0, location, length) != length) {
            throw new IllegalArgumentException("go beyond the end of the sequence");
        }
        return result;
    }
}
