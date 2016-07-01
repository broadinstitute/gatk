package org.broadinstitute.hellbender.engine.filters;

import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.utils.read.GATKRead;

/**
 * Keep only reads whose strand is forward or reverse
 */
public final class ReadStrandFilter implements ReadFilter {
    private static final long serialVersionUID = 1L;

    @Argument(fullName = "keepReverse", shortName = "kr", doc="Keep only reads on the reverse strand",optional=true)
	public boolean keepOnlyReverse = false;

    @Override
    public boolean test( final GATKRead read ) {
        return read.isReverseStrand() == keepOnlyReverse;
    }
}
