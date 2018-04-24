package org.broadinstitute.hellbender.engine.filters;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.ReadFilterArgumentDefinitions;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.Serializable;

/**
 * Keep only read pairs (0x1) with absolute insert length less than or equal to the given value.
 *
 * <p>Taking absolute values allows inclusion of pairs where the mate of the read being considered is at a smaller genomic coordinate.
 * Insert length is the difference between the 5' outer ends of mates, akin to a SAM record's TLEN (column 9).
 * Length is zero for single-end reads or when the information is unavailable.
 */
@DocumentedFeature(groupName= HelpConstants.DOC_CAT_READFILTERS, groupSummary=HelpConstants.DOC_CAT_READFILTERS_SUMMARY, summary = "Keep only read pairs with insert length less than or equal to the given value")
public final class FragmentLengthReadFilter extends ReadFilter implements Serializable  {
    private static final long serialVersionUID = 1l;

    @Argument(fullName = ReadFilterArgumentDefinitions.MAX_FRAGMENT_LENGTH_NAME,
            doc = "Maximum length of fragment (insert size)",
            optional = true)
    public int maxFragmentLength = 1000000;

    @Override
    public boolean test( final GATKRead read ) {
        if ( ! read.isPaired() ) {
            return true;
        }
        //Note fragment length is negative if mate maps to lower position than read so we take absolute value.
        return Math.abs(read.getFragmentLength()) <= maxFragmentLength;
    }
}
