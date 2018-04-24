package org.broadinstitute.hellbender.utils.smithwaterman;

import org.broadinstitute.gatk.nativebindings.smithwaterman.SWOverhangStrategy;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWParameters;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import java.io.Closeable;
import java.util.function.Supplier;

/**
 * Interface and factory for Smith-Waterman aligners
 */
public interface SmithWatermanAligner extends Closeable {
    Logger logger = LogManager.getLogger(SmithWatermanAligner.class);

    // match=1, mismatch = -1/3, gap=-(1+k/3)
    SWParameters ORIGINAL_DEFAULT = new SWParameters(3, -1, -4, -3);
    SWParameters STANDARD_NGS = new SWParameters(25, -50, -110, -6);

    /**
     *  perform a Smith-Waterman alignment of alt against ref
     *
     * @param ref bases to align to, values must be the byte equivalent of uppercase chars
     * @param alt bases to align against ref, values must be the byte equivalent of uppercase chars
     * @param parameters a set of weights to use when performing the alignment
     * @param overhangStrategy how to treat overhangs during alignment
     */
    SmithWatermanAlignment align(final byte[] ref, final byte[] alt, SWParameters parameters, SWOverhangStrategy overhangStrategy);

    /**
     * Implementations may optionally implement close in order to release any resources that they are holding.
     *
     * Calling {@link #align(byte[], byte[], SWParameters, SWOverhangStrategy)} after close must not return an incorrect
     * or invalid alignment, but otherwise the behavior is undefined.
     *
     * If a subclass implements close it's recommended that subsequent calls to align after a call to close should
     * throw {@link IllegalStateException}
     */
    @Override
    default void close() {}

    enum Implementation {
        /**
         * use the fastest available Smith-Waterman aligner that runs on your hardware
         */

        FASTEST_AVAILABLE( () -> {
            try {
                final SmithWatermanIntelAligner aligner = new SmithWatermanIntelAligner();
                logger.info("Using AVX accelerated SmithWaterman implementation");
                return aligner;
            } catch (UserException.HardwareFeatureException exception) {
                logger.info("AVX accelerated SmithWaterman implementation is not supported, falling back to the Java implementation");
                return SmithWatermanJavaAligner.getInstance();
            }
        }),

        /**
         * use the AVX enabled Smith-Waterman aligner
         */
        AVX_ENABLED( () -> {
            final SmithWatermanIntelAligner aligner = new SmithWatermanIntelAligner();
            logger.info("Using AVX accelerated SmithWaterman implementation");
            return aligner;
        }
        ),

        /**
         * use the pure java implementation of Smith-Waterman, works on all hardware
         */
        JAVA(SmithWatermanJavaAligner::getInstance);

        private final Supplier<SmithWatermanAligner> alignerSupplier;

        Implementation(final Supplier<SmithWatermanAligner> alignerSupplier ){
                this.alignerSupplier = alignerSupplier;
        }

        private SmithWatermanAligner createAligner(){
            return alignerSupplier.get();
        }
    }

    /**
     * Factory method to get an instance of an aligner corresponding to the given implementation
     */
    static SmithWatermanAligner getAligner(final Implementation type) {
        return type.createAligner();
    }
}
