package org.broadinstitute.hellbender.engine.filters;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.ReadFilterArgumentDefinitions;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import org.broadinstitute.hellbender.utils.read.GATKRead;

/**
 * Filter out reads where the number of bases without soft-clips (M, I, X, and = CIGAR operators) is lower than a threshold.
 *
 * <p>This filter is intended to filter out reads that are potentially from foreign organisms.
 * From experience with sequencing of human DNA we have found cases of contamination by bacterial
 * organisms; the symptoms of such contamination are a class of reads with only a small number
 * of aligned bases and additionally many soft-clipped bases. This filter is intended
 * to remove such reads.</p>
 *
 * <p>Note: Consecutive soft-clipped blocks are treated as a single block. For example, 1S2S10M1S2S is treated as 3S10M3S</p>
 */
@DocumentedFeature(groupName= HelpConstants.DOC_CAT_READFILTERS, groupSummary=HelpConstants.DOC_CAT_READFILTERS_SUMMARY, summary = "Filter out reads that are over-soft-clipped")
public final class OverclippedReadFilter extends ReadFilter{

    static final long serialVersionUID = 1L;

    @Argument(fullName = ReadFilterArgumentDefinitions.FILTER_TOO_SHORT_NAME,
            doc = "Minimum number of aligned bases",
            optional = true)
    public int minimumSequenceLength = 30;

    @Argument(fullName = ReadFilterArgumentDefinitions.DONT_REQUIRE_SOFT_CLIPS_BOTH_ENDS_NAME,
            doc = "Allow a read to be filtered out based on having only 1 soft-clipped block. By default, both ends must " +
                    "have a soft-clipped block, setting this flag requires only 1 soft-clipped block",
            optional = true)
    public boolean doNotRequireSoftClipsOnBothEnds;

    // Command line parser requires a no-arg constructor
    public OverclippedReadFilter() {}

    public OverclippedReadFilter(final int minimumSequenceLength, final boolean doNotRequireSoftClipsOnBothEnds) {
        this.minimumSequenceLength = minimumSequenceLength;
        this.doNotRequireSoftClipsOnBothEnds = doNotRequireSoftClipsOnBothEnds;
    }

    @Override
    public boolean test(final GATKRead read) {
        int alignedLength = 0;
        int softClipBlocks = 0;
        int minSoftClipBlocks = doNotRequireSoftClipsOnBothEnds ? 1 : 2;
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
