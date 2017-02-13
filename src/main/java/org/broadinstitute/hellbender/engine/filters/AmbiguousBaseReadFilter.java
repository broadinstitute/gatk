package org.broadinstitute.hellbender.engine.filters;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import org.broadinstitute.hellbender.utils.read.GATKRead;

/**
 * Filters out reads with more than a threshold number of N's
 */
@DocumentedFeature(groupName= HelpConstants.DOC_CAT_READFILTERS, groupSummary=HelpConstants.DOC_CAT_READFILTERS_SUMMARY)
public final class AmbiguousBaseReadFilter extends ReadFilter {

    private static final long serialVersionUID = 1L;

    @Argument(doc="Threshold fraction of non-regular bases (e.g. N) above which to filter",
            fullName="ambigFilterFrac",
            shortName="ambigFilterFrac", optional=true)
    public float N_FRAC = 0.05f;

    public AmbiguousBaseReadFilter() { }

    public AmbiguousBaseReadFilter( final float n_frac ) { this.N_FRAC = n_frac; }

    @Override
    public boolean test( final GATKRead read ) {
        final int N_max = (int)(read.getLength()*N_FRAC);
        int num_N = 0;
        for (final byte base : read.getBases()) {
            if (!BaseUtils.isRegularBase(base)) {
                num_N++;
                if (num_N > N_max) {return false;}
            }
        }
        return num_N <= N_max;
    }
}
