package org.broadinstitute.hellbender.transformers;

import htsjdk.samtools.util.SequenceUtil;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.List;

/**
 * Trims adapter sequences from read ends. If the adapter sequence ends before the end of the read, trims the remaining
 * bases as well.
 */
public class IlluminaAdapterTrimTransformer implements ReadTransformer {
    public static final long serialVersionUID = 1L;

    private final int maxMismatches;
    private final int minClipLength;
    private final List<String> adapterSequences;

    public IlluminaAdapterTrimTransformer(final int maxMismatches, final int minClipLength, final List<String> adapterSequences) {
        this.maxMismatches = maxMismatches;
        this.minClipLength = minClipLength;
        this.adapterSequences = adapterSequences;
    }

    @Override
    public GATKRead apply(final GATKRead read) {
        if (read.getLength() < minClipLength) return read;
        final byte[] bases = read.getBases();

        //Adapter list loop
        for (final String adapterSequence : adapterSequences) {
            //Start at the end of the read and walk backwards
            int alignmentLength = minClipLength;
            for (int alignmentStart = bases.length - minClipLength; alignmentStart >= 0; alignmentStart--) {
                int numMismatches = 0;
                //Check each base for mismatches
                for (int j = 0; j < alignmentLength; j++) {
                    if (!SequenceUtil.isNoCall((byte)adapterSequence.charAt(j)) && bases[alignmentStart + j] != adapterSequence.charAt(j)) {
                        if (++numMismatches > maxMismatches) break;
                    }
                }
                if (numMismatches <= maxMismatches) {
                    //We have a match. Clip from the beginning of the adapter to the end of the read.
                    final String basesString = read.getBasesString();
                    final String qualString = new String(read.getBaseQualities());
                    read.setBases(basesString.substring(0, alignmentStart).getBytes());
                    read.setBaseQualities(qualString.substring(0, alignmentStart).getBytes());
                    return read;
                }
                alignmentLength = alignmentLength < adapterSequence.length() ? alignmentLength + 1 : alignmentLength;
            }
        }
        return read;
    }

}