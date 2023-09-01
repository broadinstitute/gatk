package org.broadinstitute.hellbender.engine.filters.flow;

import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

/**
 * A common base class for flow based filters which test for conditions on an hmer basis
 */
public class FlowBasedHmerBasedReadFilterHelper  {

    interface FilterImpl {

        // provide the area of values associated with read hmers
        byte[] getValuesOfInterest(final GATKRead read);

        // check that the range of values associated with a single hmer are passing the filter
        boolean testHmer(final byte[] values, final int hmerStartingOffset, final int hmerLength);
    }

    // check if an area is a palindrome
    static boolean isPalindrome(final byte[] values, final int ofs, final int length) {

        // check that a range of bytes in the array forms an palindrome
        for (int i = 0; i < length / 2; i++) {
            if (values[ofs + i] != values[ofs + length - 1 - i]) {
                return false;
            }
        }

        return true;
    }

    static boolean test(final GATKRead read, FilterImpl impl) {

        // access qualities
        final byte[]        values = impl.getValuesOfInterest(read);
        if ( values == null )
            return false;

        // establish if edges are hard clipped
        final boolean       startHardClipped = read.getCigar().getFirstCigarElement().getOperator() == CigarOperator.HARD_CLIP;
        final boolean       endHardClipped = read.getCigar().getLastCigarElement().getOperator() == CigarOperator.HARD_CLIP;

        // iterate over hmers
        final BaseUtils.HmerIterator      iter = new BaseUtils.HmerIterator(read.getBasesNoCopy());
        int     ofs = 0;
        while ( iter.hasNext() ) {

            // find hmer
            final Pair<Byte,Integer> hmer = iter.next();
            final int                 hmerLength = hmer.getRight();

            // establish first/last
            final boolean             first = ofs == 0;
            final boolean             last = !iter.hasNext();

            // skip edge hmers if hard clipped
            if ( !((first && startHardClipped) || (last && endHardClipped)) ) {
                if (!impl.testHmer(values, ofs, hmerLength)) {
                    return false;
                }
            }

            // advance
            ofs += hmerLength;
        }

        // if here, all symetric
        return true;
    }
}
