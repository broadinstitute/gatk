package org.broadinstitute.hellbender.tools.picard.vcf.filter;

import htsjdk.samtools.util.CollectionUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFilterHeaderLine;

import java.util.*;

/**
 * Filters out sites that have a QD annotation applied to them and where the QD value is lower than a
 * lower limit.
 *
 * @author Tim Fennell
 */
public final class QdFilter implements VariantFilter {
    public static final String FILTER_NAME = "LowQD";
    private final double minimumQd;

    public QdFilter(final double minimumQd) {
        this.minimumQd = minimumQd;
    }

    @Override
    public String filter(final VariantContext ctx) {
        final double qd = ctx.getAttributeAsDouble("QD", -1d); // If QD is missing, return -1

        // QD should always be positive so a value < 0 indicates a missing value
        if (qd >= 0 && qd < minimumQd) {
            return FILTER_NAME;
        }
        else {
            return null;
        }
    }

    @Override
    public List<VCFFilterHeaderLine> headerLines() {
        return CollectionUtil.makeList(new VCFFilterHeaderLine(FILTER_NAME, "Site exhibits QD value below a hard limit."));
    }
}
