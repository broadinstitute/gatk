package org.broadinstitute.hellbender.engine.filters;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.ReadFilterArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import org.broadinstitute.hellbender.utils.read.GATKRead;

/**
 * Filter out reads where the ratio of soft-clipped bases to total bases exceeds some given value.
 *
 * Can choose to filter by number of soft-clipped bases, or by the leading/trailing soft-clipped bases only.
 */
@DocumentedFeature(
        groupName=HelpConstants.DOC_CAT_READFILTERS,
        groupSummary=HelpConstants.DOC_CAT_READFILTERS_SUMMARY,
        summary = "Filter out reads that are over-soft-clipped")
public final class SoftClippedReadFilter extends ReadFilter {
    static final long serialVersionUID = 1L;
    private final Logger logger = LogManager.getLogger(this.getClass());

    @VisibleForTesting
    @Argument(fullName = ReadFilterArgumentDefinitions.INVERT_SOFT_CLIP_RATIO_FILTER,
            doc = "Inverts the results from this filter, causing all variants that would pass to fail and visa-versa.",
            optional = true
    )
    boolean doInvertFilter = false;

    @VisibleForTesting
    @Argument(fullName = ReadFilterArgumentDefinitions.SOFT_CLIPPED_RATIO_THRESHOLD,
            doc = "Threshold ratio of soft clipped bases (anywhere in the cigar string) to total bases in read for read to be filtered.",
            optional = true,
            mutex = { ReadFilterArgumentDefinitions.SOFT_CLIPPED_LEADING_TRAILING_RATIO_THRESHOLD }
    )
    Double minimumSoftClippedRatio = null;

    @VisibleForTesting
    @Argument(fullName = ReadFilterArgumentDefinitions.SOFT_CLIPPED_LEADING_TRAILING_RATIO_THRESHOLD,
            doc = "Threshold ratio of soft clipped bases (leading / trailing the cigar string) to total bases in read for read to be filtered.",
            optional = true,
            mutex = {ReadFilterArgumentDefinitions.SOFT_CLIPPED_RATIO_THRESHOLD}
    )
    Double minimumLeadingTrailingSoftClippedRatio = null;

    // Command line parser requires a no-arg constructor
    public SoftClippedReadFilter() {}

    private boolean testMinSoftClippedRatio(final GATKRead read) {
        int totalLength = 0;
        int numSoftClippedBases = 0;
        for ( final CigarElement element : read.getCigarElements() ) {
            if (element.getOperator() == CigarOperator.S) {
                numSoftClippedBases += element.getLength();
            }
            totalLength += element.getLength();
        }

        final double softClipRatio = ((double)numSoftClippedBases / (double)totalLength);

        return softClipRatio > minimumSoftClippedRatio;
    }

    private boolean testMinLeadingTrailingSoftClippedRatio(final GATKRead read) {

        if ( read.getCigarElements().size() < 1 ) {
            return false;
        }

        // Get the index of the last cigar element:
        final int lastCigarOpIndex = read.getCigarElements().size() - 1;

        // Calculate the number of bases here:
        // Note: we don't want to double-count if the read has only 1 cigar operator.
        final int numLeadingTrailingSoftClippedBases =
                (read.getCigarElement(0).getOperator() == CigarOperator.S ? read.getCigarElement(0).getLength() : 0)
                        +
                ((lastCigarOpIndex != 0) && (read.getCigarElement(lastCigarOpIndex).getOperator() == CigarOperator.S)
                        ? read.getCigarElement(lastCigarOpIndex).getLength()
                        : 0);

        // Determine length:
        final int totalLength = read.getCigarElements().stream()
                .mapToInt( CigarElement::getLength )
                .sum();

        // Calculate the ratio:
        final double softClipRatio = ((double)numLeadingTrailingSoftClippedBases / (double)totalLength);

        return softClipRatio > minimumLeadingTrailingSoftClippedRatio;
    }

    @Override
    public boolean test(final GATKRead read) {

        final boolean result;

        // NOTE: Since we have mutex'd the args for the clipping ratios, we only need to see if they
        //       have been specified.  If they have, that's the filter logic we're using.
        // If we specified the clipping ratio, we use the min sequence length test:
        if ( minimumSoftClippedRatio != null ) {
            result = testMinSoftClippedRatio(read);
        }
        // If we specified the leading/trailing clipping ratio, we use the min sequence length test:
        else if ( minimumLeadingTrailingSoftClippedRatio != null ) {
            result = testMinLeadingTrailingSoftClippedRatio(read);
        }
        else {
            throw new UserException(
                    "Must provide one of the following arguments: " +
                            ReadFilterArgumentDefinitions.SOFT_CLIPPED_RATIO_THRESHOLD + "," +
                            ReadFilterArgumentDefinitions.SOFT_CLIPPED_LEADING_TRAILING_RATIO_THRESHOLD
            );
        }

        // Check for if we want to invert our results:
        if ( doInvertFilter ) {
            return !result;
        }
        return result;
    }
}
