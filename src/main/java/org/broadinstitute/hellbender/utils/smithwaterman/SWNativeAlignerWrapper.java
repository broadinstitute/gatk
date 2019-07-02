package org.broadinstitute.hellbender.utils.smithwaterman;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWAlignerNativeBinding;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWNativeAlignerResult;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWOverhangStrategy;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWParameters;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.AlignmentUtils;

import java.util.ArrayList;
import java.util.List;

/**
 * A wrapper that converts instances of {@link SWAlignerNativeBinding} into a {@link SmithWatermanAligner}
 */
public final class SWNativeAlignerWrapper implements SmithWatermanAligner {
    private final SWAlignerNativeBinding aligner;
    private long totalComputeTime = 0;

    public SWNativeAlignerWrapper(final SWAlignerNativeBinding aligner) {
        this.aligner = aligner;
    }

    @Override
    public SmithWatermanAlignment align(final byte[] reference, final byte[] alternate, final SWParameters parameters, final SWOverhangStrategy overhangStrategy){
        long startTime = System.nanoTime();

        Utils.nonNull(parameters);
        Utils.nonNull(overhangStrategy);

        // avoid running full Smith-Waterman if there is an exact match of alternate in reference
        int matchIndex = -1;
        if (overhangStrategy == SWOverhangStrategy.SOFTCLIP || overhangStrategy == SWOverhangStrategy.IGNORE) {
            // Use a substring search to find an exact match of the alternate in the reference
            // NOTE: This approach only works for SOFTCLIP and IGNORE overhang strategies
            matchIndex = Utils.lastIndexOf(reference, alternate);
        }

        final SmithWatermanAlignment alignmentResult;

        if (matchIndex != -1) {
            // generate the alignment result when the substring search was successful
            final List<CigarElement> lce = new ArrayList<>(alternate.length);
            lce.add(new CigarElement(alternate.length, CigarOperator.M));
            alignmentResult = new SWNativeResultWrapper(AlignmentUtils.consolidateCigar(new Cigar(lce)), matchIndex);
        }
        else {
            // run full Smith-Waterman
            final SWNativeAlignerResult alignment = aligner.align(reference, alternate,parameters,overhangStrategy);
            alignmentResult =  new SWNativeResultWrapper(alignment);
        }

        totalComputeTime += System.nanoTime() - startTime;
        return alignmentResult;
    }

    /**
     * Report total compute time and close aligner
     */
    @Override
    public void close() {
        logger.info(String.format("Total compute time in native Smith-Waterman : %.2f sec", totalComputeTime * 1e-9));
        aligner.close();
    }

    private static final class SWNativeResultWrapper implements SmithWatermanAlignment {
        private final Cigar cigar;
        private final int alignmentOffset;

        public SWNativeResultWrapper(final SWNativeAlignerResult nativeResult) {
            this.cigar = TextCigarCodec.decode(nativeResult.cigar);
            this.alignmentOffset = nativeResult.alignment_offset;
        }

        public SWNativeResultWrapper(final Cigar cigar, final int alignmentOffset) {
            this.cigar = cigar;
            this.alignmentOffset =  alignmentOffset;
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
