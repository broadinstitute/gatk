package org.broadinstitute.hellbender.utils.smithwaterman;

import org.broadinstitute.gatk.nativebindings.smithwaterman.SWOverhangStrategy;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWParameters;

import java.io.Closeable;
import java.util.function.Supplier;

/**
 * Interface for Smith-Waterman aligners
 */
public interface SmithWatermanAligner extends Closeable {

    // match=1, mismatch = -1/3, gap=-(1+k/3)
    SWParameters ORIGINAL_DEFAULT = new SWParameters(3, -1, -4, -3);
    SWParameters STANDARD_NGS = new SWParameters(25, -50, -110, -6);

    /**
      *  Perform an alignment and return the result
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
        FASTEST_AVAILABLE(SmithWatermanJavaAligner::getInstance),
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
