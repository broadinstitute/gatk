package org.broadinstitute.hellbender.utils.smithwaterman;

import org.broadinstitute.gatk.nativebindings.smithwaterman.SWOverhangStrategy;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWParameters;

import java.io.Closeable;
import java.util.function.Supplier;

public interface SmithWatermanAligner extends Closeable{

    SmithWatermanAlignment align(final byte[] ref, final byte[] alt, SWParameters parameters, SWOverhangStrategy overhangStrategy);

    @Override
    default void close() {}

    enum Implementation {
        FASTEST_AVAILABLE(SmithWatermanJavaAligner::getInstance),
        JAVA(SmithWatermanJavaAligner::getInstance);

        private final Supplier<SmithWatermanAligner> createAligner;

        Implementation(final Supplier<SmithWatermanAligner> createAligner ){
                this.createAligner = createAligner;
        }

        private SmithWatermanAligner createAligner(){
            return createAligner.get();
        }
    }

    static SmithWatermanAligner getAligner(final Implementation type) {
        return type.createAligner();
    }
}
