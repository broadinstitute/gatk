package org.broadinstitute.hellbender.tools.spark.sv.utils;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.CigarBuilder;
import org.broadinstitute.hellbender.utils.read.CigarUtils;
import org.broadinstitute.hellbender.utils.read.ClippingTail;

import java.util.ArrayList;
import java.util.List;

/**
 * Various utility functions helping calling structural variants.
 */
public final class SvCigarUtils {

    /**
     * Convert the 'I' CigarElement, if it is at either end (terminal) of the input cigar, to a corresponding 'S' operator.
     * Note that we allow CIGAR of the format '10H10S10I10M', but disallows the format if after the conversion the cigar turns into a giant clip,
     * e.g. '10H10S10I10S10H' is not allowed (if allowed, it becomes a giant clip of '10H30S10H' which is non-sense).
     *
     * @return a pair of number of clipped (hard and soft, including the ones from the converted terminal 'I') bases at the front and back of the
     *         input {@code cigarAlongInput5to3Direction}.
     *
     * @throws IllegalArgumentException when the checks as described above fail.
     */
    public static Cigar convertTerminalInsertionToSoftClip(final Cigar cigar) {

        if (cigar.numCigarElements() < 2 ) {
            return cigar;
        }

        final CigarBuilder builder = new CigarBuilder();
        for (int n = 0; n < cigar.numCigarElements(); n++) {
            final CigarElement element = cigar.getCigarElement(n);
            if (element.getOperator() != CigarOperator.INSERTION) { // not an insertion
                builder.add(element);
            } else if (n == 0 || n == cigar.numCigarElements() - 1) {   // terminal insertion with no clipping -- convert to soft clip
                builder.add(new CigarElement(element.getLength(), CigarOperator.SOFT_CLIP));
            } else if (cigar.getCigarElement(n-1).getOperator().isClipping() || cigar.getCigarElement(n+1).getOperator().isClipping()) {    // insertion preceding or following clip
                builder.add(new CigarElement(element.getLength(), CigarOperator.SOFT_CLIP));
            } else {    // interior insertion
                builder.add(element);
            }
        }

        return builder.make();
    }


}
