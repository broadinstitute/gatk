package org.broadinstitute.hellbender.engine.filters;

import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.cmdline.Argument;

/**
 * Keep only reads whose strand is negative or positive
 */
public final class ReadStrandFilter implements ReadFilter {

    @Argument(fullName = "keepReverse", shortName = "kr", doc="Keep only reads on the reverse strand",optional=true)
	boolean keepOnlyReverse = false;

    @Override
    public boolean test(final SAMRecord read) {
        return read.getReadNegativeStrandFlag() == keepOnlyReverse;
    }
}
