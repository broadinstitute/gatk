package org.broadinstitute.hellbender.engine.filters;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import org.broadinstitute.hellbender.utils.read.GATKRead;

/**
 * Filters out reads with more than a threshold number of N's. By default will determine the threshold by the given
 * maximum fraction of bases in the read. This can be overridden by setting the maximum number of ambiguous bases instead.
 */
@DocumentedFeature(groupName= HelpConstants.DOC_CAT_READFILTERS, groupSummary=HelpConstants.DOC_CAT_READFILTERS_SUMMARY)
public final class AmbiguousBaseReadFilter extends ReadFilter {

    private static final long serialVersionUID = 1L;
    private final Integer maxAmbiguousBases;

    @Argument(doc = "Threshold fraction of ambiguous bases",
            fullName = "ambigFilterFrac",
            optional = true)
    public double maxAmbiguousBaseFraction = 0.05;

    public AmbiguousBaseReadFilter() { maxAmbiguousBases = null; }

    /**
     * Sets max number of ambiguous bases, which overrides maxAmbiguousBaseFraction
     */
    public AmbiguousBaseReadFilter(final int maxAmbiguousBases) {
        this.maxAmbiguousBases = maxAmbiguousBases;
    }

    //Filters out reads with more than a threshold number of N's
    @Override
    public boolean test(final GATKRead read) {
        final int maxN = maxAmbiguousBases != null ? maxAmbiguousBases : (int) (read.getLength() * maxAmbiguousBaseFraction);
        int numN = 0;
        for (final byte base : read.getBases()) {
            if (!BaseUtils.isRegularBase(base)) {
                numN++;
                if (numN > maxN) {
                    return false;
                }
            }
        }
        return true;
    }
}
