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

    @Argument(doc = "Threshold fraction of ambiguous bases",
            fullName = "ambigFilterFrac",
            optional = true)
    public double nFrac = 0.05;

    public AmbiguousBaseReadFilter() {
    }

    public AmbiguousBaseReadFilter(final double n_frac) {
        this.nFrac = n_frac;
    }

    //Filters out reads with more than a threshold number of N's
    @Override
    public boolean test(final GATKRead read) {
        final int N_max = (int) (read.getLength() * nFrac);
        int num_N = 0;
        for (final byte base : read.getBases()) {
            if (!BaseUtils.isRegularBase(base)) {
                num_N++;
                if (num_N > N_max) {
                    return false;
                }
            }
        }
        return num_N <= N_max;
    }
}
