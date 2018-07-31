package org.broadinstitute.hellbender.engine.filters;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.ReadFilterArgumentDefinitions;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import org.broadinstitute.hellbender.utils.read.GATKRead;

/**
 * Keep only reads whose strand is either forward (not 0x10) or reverse (0x10), as specified.
 *
 * <p>By default the filter keeps only forward reads (not 0x10).</p>
 */
@DocumentedFeature(groupName= HelpConstants.DOC_CAT_READFILTERS, groupSummary=HelpConstants.DOC_CAT_READFILTERS_SUMMARY, summary = "Keep only reads whose strand is as specified")
public final class ReadStrandFilter extends ReadFilter {
    private static final long serialVersionUID = 1L;

    @Argument(fullName = ReadFilterArgumentDefinitions.KEEP_REVERSE_STRAND_ONLY_NAME,
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
