package org.broadinstitute.hellbender.transformers;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.concurrent.atomic.AtomicInteger;

/**
 * Checks for and errors out (or fixes if requested) when it detects reads with base qualities that are not encoded with
 * phred-scaled quality scores.  Q0 == ASCII 33 according to the SAM specification, whereas Illumina encoding starts at
 * Q64.  The idea here is simple: if we are asked to fix the scores then we just subtract 31 from every quality score.
 */
public final class MisencodedBaseQualityReadTransformer implements ReadTransformer {
    private static final long serialVersionUID = 1L;

    private static final int ILLUMINA_ENCODING_FIX_VALUE = 31;  // Illumina_64 - PHRED_33

    private final boolean fixQuals;

    protected static final int samplingFrequency = 1000; // sample 1 read for every 1000 encountered
    // atomic integer in case of concurrent usage of this transformer
    protected static AtomicInteger currentReadCounter = new AtomicInteger(0);

    /**
     * Constructor for allow only checking of misencoded qualities
     *
     * @param fixQuals if {@code false} it will just check the qualities; otherwise, it will fix them
     */
    public MisencodedBaseQualityReadTransformer(final boolean fixQuals) {
        this.fixQuals = fixQuals;
    }

    /**
     * Constructor for fixing misencoded qualities
     */
    public MisencodedBaseQualityReadTransformer() {
        this(true);
    }

    @Override
    public GATKRead apply(final GATKRead read) {
        if(fixQuals) {
            fixMisencodedQuals(read);
        } else {
            checkForMisencodedQuals(read);
        }
        return read;
    }

    private void fixMisencodedQuals(final GATKRead read) {
        final byte[] quals = read.getBaseQualities();
        for ( int i = 0; i < quals.length; i++ ) {
            quals[i] -= ILLUMINA_ENCODING_FIX_VALUE;
            if ( quals[i] < 0 )
                throw new UserException.BadInput("while fixing mis-encoded base qualities we encountered a read that was correctly encoded; we cannot handle such a mixture of reads so unfortunately the BAM must be fixed with some other tool");
        }
        read.setBaseQualities(quals);
    }

    private void checkForMisencodedQuals(final GATKRead read) {
        // sample reads randomly for checking
        if ( currentReadCounter.incrementAndGet() >= samplingFrequency ) {
            currentReadCounter.set(0);
            final byte[] quals = read.getBaseQualities();
            for ( final byte qual : quals ) {
                if ( qual > QualityUtils.MAX_REASONABLE_Q_SCORE )
                    throw new UserException.MisencodedQualityScoresRead(read, "we encountered an extremely high quality score of " + (int) qual);
            }
        }
    }

}
