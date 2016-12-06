package org.broadinstitute.hellbender.engine.filters;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.utils.read.GATKRead;

/**
 * Keep only reads whose strand is forward or reverse
 */
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
