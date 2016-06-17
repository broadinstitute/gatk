package org.broadinstitute.hellbender.engine.filters;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.cmdline.Argument;

/**
 * Filter out reads that are over-soft-clipped
 *
 * <p>
 *     This filter is intended to filter out reads that are potentially from foreign organisms.
 *     From experience with sequencing of human DNA we have found cases of contamination by bacterial
 *     organisms; the symptoms of such contamination are a class of reads with only a small number
 *     of aligned bases and additionally many soft-clipped bases.  This filter is intended
 *     to remove such reads. Consecutive soft-clipped blocks are treated as a single block
 * </p>
 *
 */
public final class OverclippedReadFilter implements ReadFilter{

    static final long serialVersionUID = 1L;

    @Argument(fullName = "filter_is_too_short_value", shortName = "filterTooShort", doc = "Value for which reads with less than this number of aligned bases is considered too short", optional = true)
    int minimumSequenceLength = 30;

    @Argument(fullName = "do_not_require_softclips_both_ends", shortName = "NoRequireSCBothEnds", doc = "Allow a read to be filtered out based on having only 1 soft-clipped block. By default, both ends must have a soft-clipped block, setting this flag requires only 1 soft-clipped block.", optional = true)
    Boolean doNotRequireSoftclipsOnBothEnds = false;


    @Override
    public boolean test(final GATKRead read) {
        int alignedLength = 0;
        int softClipBlocks = 0;
        int minSoftClipBlocks = doNotRequireSoftclipsOnBothEnds ? 1 : 2;
        CigarOperator prevOperator = null;

        for ( final CigarElement element : read.getCigarElements() ) {
            if ( element.getOperator() == CigarOperator.S ) { // count total S blocks
                //Treat consecutive S blocks as a single one
                if(prevOperator != CigarOperator.S){
                    softClipBlocks += 1;
                }

            } else if ( element.getOperator().consumesReadBases() ) {   // M, I, X, and EQ (S was already accounted for above)
                alignedLength += element.getLength();
            }
            prevOperator = element.getOperator();
        }

        return(alignedLength >= minimumSequenceLength || softClipBlocks < minSoftClipBlocks);
    }
}
