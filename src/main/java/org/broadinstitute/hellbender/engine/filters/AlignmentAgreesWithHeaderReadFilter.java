package org.broadinstitute.hellbender.engine.filters;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.io.Serializable;

/**
 * Filter out reads where the alignment does not match the contents of the header.
 *
 * <p>The read does not match the contents of the header if:</p>
 *
 * <ul>
 *     <li>It is aligned to a non-existing contig</li>
 *     <li>It is aligned to a point after the end of the contig</li>
 * </ul>
 */
@DocumentedFeature(groupName= HelpConstants.DOC_CAT_READFILTERS, groupSummary=HelpConstants.DOC_CAT_READFILTERS_SUMMARY, summary = "Filters out reads where the alignment does not match the contents of the header")
public final class AlignmentAgreesWithHeaderReadFilter extends ReadFilter implements Serializable {
    private static final long serialVersionUID = 1l;

    public AlignmentAgreesWithHeaderReadFilter( ) {}

    public AlignmentAgreesWithHeaderReadFilter( final SAMFileHeader header ) {
        super.setHeader(header);
    }

    @Override
    public boolean test( GATKRead read ) {
        return ReadUtils.alignmentAgreesWithHeader(samHeader, read);
    }
}
