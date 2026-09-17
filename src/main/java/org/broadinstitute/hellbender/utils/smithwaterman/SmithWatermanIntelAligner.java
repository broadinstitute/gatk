package org.broadinstitute.hellbender.utils.smithwaterman;

import com.fulcrumgenomics.fgkl.smithwaterman.FgklSmithWaterman;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWOverhangStrategy;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWParameters;
import org.broadinstitute.hellbender.exceptions.UserException;

/**
 * A {@link SmithWatermanAligner} backed by the fgkl native library, which selects the fastest SIMD kernel available
 * on the running CPU. Construction throws {@link UserException.HardwareFeatureException} when no native library is
 * available for this platform.
 */
public final class SmithWatermanIntelAligner implements SmithWatermanAligner {

    private static final Logger logger = LogManager.getLogger(SmithWatermanIntelAligner.class);

    private final FgklSmithWaterman smithWaterman = new FgklSmithWaterman();
    private final SWNativeAlignerWrapper alignerWrapper = new SWNativeAlignerWrapper(smithWaterman);

    public SmithWatermanIntelAligner() throws UserException.HardwareFeatureException {
        if (!smithWaterman.load(null)) {
            throw new UserException.HardwareFeatureException("The native Smith-Waterman library is not available on this platform.");
        }
        logger.info("Using the native Smith-Waterman aligner with the " + smithWaterman.backend() + " backend");
    }

    /**
     * Aligns the alternate sequence to the reference sequence
     *
     * @param reference  ref sequence
     * @param alternate  alt sequence
     */
    @Override
    public SmithWatermanAlignment align(final byte[] reference, final byte[] alternate, final SWParameters parameters, final SWOverhangStrategy overhangStrategy) {
        return alignerWrapper.align(reference, alternate, parameters, overhangStrategy);
    }

    /**
     * Close the aligner
     */
    @Override
    public void close() {
        alignerWrapper.close();
    }
}
