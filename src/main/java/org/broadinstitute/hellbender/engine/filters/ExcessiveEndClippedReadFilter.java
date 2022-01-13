package org.broadinstitute.hellbender.engine.filters;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.ReadFilterArgumentDefinitions;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.Collections;
import java.util.List;

/**
 * Filter out reads where the number of soft-/hard-clipped bases on either end is above a certain threshold.
 *
 * This filter will filter reads where the number of soft clips on one end of the read is too large.  It does NOT
 * check the total number of bases clipped on BOTH ends.  For the purposes of this filter, the count of soft and hard
 * clipped bases on either end of a read are combined.
 *
 * For example, a read with the following cigar strings will be filtered out (given maxClippedBases=1000):
 *
 *      1500S0000M
 *      900H300S50000M
 *      1001H50000M
 *
 * However, this read will NOT be filtered out (given maxClippedBases=1000):
 *
 *      900S50000M900S
 *      150H800S50000M900S150H
 *      1000S500M
 *      500M1000S
 *
 * Note: This is designed for use with (unpaired) long reads, but can be used with any reads.
 */
@DocumentedFeature(
        groupName= HelpConstants.DOC_CAT_READFILTERS,
        groupSummary=HelpConstants.DOC_CAT_READFILTERS_SUMMARY,
        summary = "Filter out reads that have too many clipped bases on either end."
)
public final class ExcessiveEndClippedReadFilter extends ReadFilter {
    static final long serialVersionUID = 1L;
    private final Logger logger = LogManager.getLogger(this.getClass());

    @Argument(fullName = ReadFilterArgumentDefinitions.MAXIMUM_END_CLIPPED_BASES,
            doc = "Maximum number of clipped bases on either end of a given read",
            optional = true)
    private int maxClippedBases = 1000;

    // Command line parser requires a no-arg constructor
    public ExcessiveEndClippedReadFilter() {}

    @VisibleForTesting
    ExcessiveEndClippedReadFilter(final int maxClippedBases) {
        this.maxClippedBases = maxClippedBases;
    }

    @Override
    public boolean test(final GATKRead read) {
        final List<CigarElement> cigarElements = read.getCigarElements();

        // Check front end of the read:
        int clippedBases = 0;
        for ( final CigarElement element : cigarElements) {
            if ( element.getOperator() != CigarOperator.S && element.getOperator() != CigarOperator.H) {
                break;
            }
            clippedBases += element.getLength();
        }
        // If we fail on the front end we can quit now:
        if (clippedBases > maxClippedBases) {
            return false;
        }

        // Check back end of the read:
        clippedBases = 0;
        for ( int i = cigarElements.size() - 1; i >= 0; --i ) {
            final CigarElement element = cigarElements.get(i);
            if ( element.getOperator() != CigarOperator.S && element.getOperator() != CigarOperator.H) {
                break;
            }
            clippedBases += element.getLength();
        }
        return clippedBases <= maxClippedBases;
    }
}
