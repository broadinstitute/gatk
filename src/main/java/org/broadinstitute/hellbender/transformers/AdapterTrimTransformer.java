package org.broadinstitute.hellbender.transformers;

import htsjdk.samtools.util.SequenceUtil;
import org.broadinstitute.hellbender.utils.clipping.ClippingOp;
import org.broadinstitute.hellbender.utils.clipping.ClippingRepresentation;
import org.broadinstitute.hellbender.utils.clipping.ReadClipper;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.List;

/**
 * Trims (hard clips) adapter sequences from read ends. If the adapter sequence ends before the end of the read,
 * trims the remaining bases as well. Trims the closest-matching adapter sequence. In the case of a tie, trims
 * the adapter that occurs earliest in the read.
 */
public final class AdapterTrimTransformer implements ReadTransformer {
    private static final long serialVersionUID = 1L;

    private final int maxMismatches;
    private final int minClipLength;
    private final List<String> adapterSequences;

    public AdapterTrimTransformer(final int maxMismatches, final int minClipLength, final List<String> adapterSequences) {
        this.maxMismatches = maxMismatches;
        this.minClipLength = minClipLength;
        this.adapterSequences = adapterSequences;
    }

    @Override
    public GATKRead apply(final GATKRead read) {
        if (read.getLength() < minClipLength) return read;
        final byte[] bases = read.getBases();

        //Be sure to find the best match in case one adapter ends with another adapter's sequence (within maxMismatches difference)
        int bestNumMismatches = maxMismatches;
        int bestAlignmentStart = bases.length;
        boolean foundAdapter = false;

        //Adapter list loop
        for (final String adapterSequence : adapterSequences) {
            //Start at the end of the read and walk backwards
            int alignmentLength = minClipLength;
            for (int alignmentStart = bases.length - minClipLength; alignmentStart >= 0; alignmentStart--) {
                int numMismatches = 0;
                //Check each base for mismatches
                for (int j = 0; j < alignmentLength; j++) {
                    if (!SequenceUtil.isNoCall((byte) adapterSequence.charAt(j)) && bases[alignmentStart + j] != adapterSequence.charAt(j)) {
                        if (++numMismatches > maxMismatches) break;
                    }
                }
                if (numMismatches < bestNumMismatches || (numMismatches == bestNumMismatches && alignmentStart < bestAlignmentStart)) {
                    //We have a (better/earlier) match
                    bestNumMismatches = numMismatches;
                    bestAlignmentStart = alignmentStart;
                    foundAdapter = true;
                }
                alignmentLength = alignmentLength < adapterSequence.length() ? alignmentLength + 1 : alignmentLength;
            }
        }
        if (foundAdapter) {
            //Hard clip from the beginning of the adapter to the end of the read
            final ReadClipper readClipper = new ReadClipper(read);
            readClipper.addOp(new ClippingOp(bestAlignmentStart, read.getLength()));
            return readClipper.clipRead(ClippingRepresentation.HARDCLIP_BASES);
        }
        return read;
    }

}