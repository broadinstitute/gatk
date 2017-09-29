package org.broadinstitute.hellbender.engine.filters;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import org.broadinstitute.hellbender.utils.read.GATKRead;

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
@DocumentedFeature(groupName= HelpConstants.DOC_CAT_READFILTERS, groupSummary=HelpConstants.DOC_CAT_READFILTERS_SUMMARY)
public final class OverclippedReadFilter extends ReadFilter{

    static final long serialVersionUID = 1L;

    @Argument(fullName = "filterTooShort",
            shortName = "filterTooShort",
            doc = "Value for which reads with less than this number of aligned bases is considered too short",
            optional = true)
    public int minimumSequenceLength = 30;

    @Argument(fullName = "dontRequireSoftClipsBothEnds",
            shortName = "dontRequireSoftClipsBothEnds",
            doc = "Allow a read to be filtered out based on having only 1 soft-clipped block. By default, both ends must " +
                    "have a soft-clipped block, setting this flag requires only 1 soft-clipped block.",
            optional = true)
    public boolean doNotRequireSoftclipsOnBothEnds;

    // Command line parser requires a no-arg constructor
    public OverclippedReadFilter() {}

    public OverclippedReadFilter(final int minimumSequenceLength, final boolean doNotRequireSoftclipsOnBothEnds) {
        this.minimumSequenceLength = minimumSequenceLength;
        this.doNotRequireSoftclipsOnBothEnds = doNotRequireSoftclipsOnBothEnds;
    }

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
