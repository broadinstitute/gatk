package org.broadinstitute.hellbender.engine.filters;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.ReadFilterArgumentDefinitions;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import org.broadinstitute.hellbender.utils.read.GATKRead;

/**
 * Filters out reads that have greater than the threshold number for unknown (N) bases.
 */
@DocumentedFeature(groupName= HelpConstants.DOC_CAT_READFILTERS, groupSummary=HelpConstants.DOC_CAT_READFILTERS_SUMMARY, summary = "Filters out reads that have greater than the threshold number of N bases")
public final class AmbiguousBaseReadFilter extends ReadFilter {

    private static final long serialVersionUID = 1L;

    @Argument(doc = "Threshold fraction of ambiguous bases",
            fullName = ReadFilterArgumentDefinitions.AMBIGUOUS_FILTER_FRACTION_NAME,
            optional = true,
            mutex = {ReadFilterArgumentDefinitions.AMBIGUOUS_FILTER_BASES_NAME})
    public double maxAmbiguousBaseFraction = 0.05;

    @Argument(doc = "Threshold number of ambiguous bases. If null, uses threshold fraction; otherwise, overrides threshold fraction.",
            fullName = ReadFilterArgumentDefinitions.AMBIGUOUS_FILTER_BASES_NAME,
            optional = true,
            mutex = {ReadFilterArgumentDefinitions.AMBIGUOUS_FILTER_FRACTION_NAME})
    public Integer maxAmbiguousBases = null;

    public AmbiguousBaseReadFilter() {}

    public AmbiguousBaseReadFilter(final int maxAmbiguousBases) {
        this.maxAmbiguousBases = maxAmbiguousBases;
    }

    /**
     * Test read for a given maximum threshold of allowable ambiguous bases
     */
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
