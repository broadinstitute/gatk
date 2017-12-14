package org.broadinstitute.hellbender.utils.smithwaterman;

import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWAlignerNativeBinding;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWNativeAlignerResult;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWOverhangStrategy;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWParameters;

import java.io.File;

public class SWNativeAlignerWrapperUnitTest extends SmithWatermanAlignerAbstractUnitTest {


    @Override
    protected SmithWatermanAligner getAligner() {
        final SWAlignerNativeBinding javaBackedNativeBinding = new SWAlignerNativeBinding() {

            @Override
            public boolean load(File tmpDir) {
                return true;
            }

            @Override
            public SWNativeAlignerResult align(byte[] ref, byte[] alt, SWParameters parameters, SWOverhangStrategy overhangStrategy) {
                final SmithWatermanAlignment alignment = SmithWatermanJavaAligner.getInstance().align(ref, alt,
                                                                                                          parameters,
                                                                                                          overhangStrategy);
                return new SWNativeAlignerResult(TextCigarCodec.encode(alignment.getCigar()), alignment.getAlignmentOffset());
            }
        };
        return new SWNativeAlignerWrapper(javaBackedNativeBinding);
    }
}
