package org.broadinstitute.hellbender.engine.filters;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import org.broadinstitute.hellbender.utils.read.GATKRead;

/**
 * Keep only reads whose strand is forward or reverse
 */
@DocumentedFeature(groupName= HelpConstants.DOC_CAT_READFILTERS, groupSummary=HelpConstants.DOC_CAT_READFILTERS_SUMMARY)
public final class ReadStrandFilter extends ReadFilter {
    private static final long serialVersionUID = 1L;

    @Argument(fullName = "keepReverse",
            shortName = "keepReverse",
            doc="Keep only reads on the reverse strand",
            optional=false)
	public Boolean keepOnlyReverse;

    public ReadStrandFilter() {}

    public ReadStrandFilter(final boolean keepOnlyReverse) { this.keepOnlyReverse = keepOnlyReverse; }

    @Override
    public boolean test( final GATKRead read ) {
        return read.isReverseStrand() == keepOnlyReverse;
    }
}
