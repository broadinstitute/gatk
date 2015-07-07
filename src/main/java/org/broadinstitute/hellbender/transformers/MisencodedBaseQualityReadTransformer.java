package org.broadinstitute.hellbender.transformers;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.read.GATKRead;

/**
 * Checks for and errors out (or fixes if requested) when it detects reads with base qualities that are not encoded with
 * phred-scaled quality scores.  Q0 == ASCII 33 according to the SAM specification, whereas Illumina encoding starts at
 * Q64.  The idea here is simple: if we are asked to fix the scores then we just subtract 31 from every quality score.
 */
public final class MisencodedBaseQualityReadTransformer implements ReadTransformer {
    private static final long serialVersionUID = 1L;

    private static final int ILLUMINA_ENCODING_FIX_VALUE = 31;  // Illumina_64 - PHRED_33

    @Override
    public GATKRead apply(final GATKRead read) {
        final byte[] quals = read.getBaseQualities();
        for ( int i = 0; i < quals.length; i++ ) {
            quals[i] -= ILLUMINA_ENCODING_FIX_VALUE;
            if ( quals[i] < 0 )
                throw new UserException.BadInput("while fixing mis-encoded base qualities we encountered a read that was correctly encoded; we cannot handle such a mixture of reads so unfortunately the BAM must be fixed with some other tool");
        }
        read.setBaseQualities(quals);
        return read;
    }
}
