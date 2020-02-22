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

    default byte[] bytesAt(long location, final int length) {
        final byte[] result = new byte[length];

        for (int i = 0; i < length; i++, location++) {
            result[i] = byteAt(location);
        }
        return result;
    }
}
