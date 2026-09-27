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
         * use the native Smith-Waterman aligner when a native library is available for this platform, and the Java
         * aligner otherwise
         */
        FASTEST_AVAILABLE( () -> {
            try {
                return new SmithWatermanNativeAligner();
            } catch (UserException.HardwareFeatureException exception) {
                logger.warn("***********************************************************************************************");
                logger.warn("*** WARNING: no native Smith-Waterman library is available for this platform: " + exception.getMessage());
                logger.warn("*** Falling back to the slower Java Smith-Waterman aligner!");
                logger.warn("***********************************************************************************************");
                return SmithWatermanJavaAligner.getInstance();
            }
        }),

        /**
         * use the native Smith-Waterman aligner; fails if no native library is available for this platform
         */
        AVX_ENABLED(SmithWatermanNativeAligner::new),

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
