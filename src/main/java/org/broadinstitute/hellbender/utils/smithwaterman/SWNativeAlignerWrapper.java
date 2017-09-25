package org.broadinstitute.hellbender.utils.smithwaterman;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWAlignerNativeBinding;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWNativeAlignerResult;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWOverhangStrategy;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWParameters;

/**
 * A wrapper that converts instances of {@link SWAlignerNativeBinding} into a {@link SmithWatermanAligner}
 */
public final class SWNativeAlignerWrapper implements SmithWatermanAligner {
    private final SWAlignerNativeBinding aligner;

    public SWNativeAlignerWrapper(final SWAlignerNativeBinding aligner) {
        this.aligner = aligner;
    }

    @Override
    public SmithWatermanAlignment align(final byte[] ref, final byte[] alt, final SWParameters parameters, final SWOverhangStrategy overhangStrategy){
        final SWNativeAlignerResult alignment = aligner.align(ref, alt, parameters, overhangStrategy);
        return new SWNativeResultWrapper(alignment);
    }

    @Override
    public void close() {
        aligner.close();
    }

    private static final class SWNativeResultWrapper implements SmithWatermanAlignment {
        private final Cigar cigar;
        private final int alignmentOffset;

        public SWNativeResultWrapper(final SWNativeAlignerResult nativeResult) {
            this.cigar = TextCigarCodec.decode(nativeResult.cigar);
            this.alignmentOffset = nativeResult.alignment_offset;
        }

        @Override
        public Cigar getCigar() {
            return cigar;
        }

        @Override
        public int getAlignmentOffset() {
            return alignmentOffset;
        }
    }
}
