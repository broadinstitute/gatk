package org.broadinstitute.hellbender.transformers;

import org.broadinstitute.hellbender.utils.read.GATKRead;

/**
 * Masks read bases with a supra-threshold number of A/T's or G/C's within a given window size.
 */
public class SimpleRepeatMaskTransformer implements ReadTransformer {
    public static final long serialVersionUID = 1L;
    private static final byte MASKED_BASE_QUAL = 2;

    final int windowSize, threshAT, threshGC;

    public SimpleRepeatMaskTransformer(final int threshAT, final int threshGC, final int windowSize) {
        this.threshAT = threshAT;
        this.threshGC = threshGC;
        this.windowSize = windowSize;
    }

    @Override
    public GATKRead apply(final GATKRead read) {
        final byte[] bases = read.getBases();
        final byte[] quals = read.getBaseQualities();
        final byte[] originalBases = bases.clone();
        byte[] baseCounts = new byte[255];

        //Initialize with the first window
        int windowStart = 0;
        int windowEnd = Math.min(windowSize - 1, bases.length - 1);
        for (int i = windowStart; i <= windowEnd; i++) {
            baseCounts[originalBases[i]]++;
        }
        maskWindow(baseCounts, windowStart, windowEnd, bases, quals);

        //March the window through the read, one base at a time
        while (windowEnd < bases.length - 1) {
            baseCounts[originalBases[windowStart++]]--; //Subtract old first base
            baseCounts[originalBases[++windowEnd]]++; //Add new final base
            maskWindow(baseCounts, windowStart, windowEnd, bases, quals);
        }
        read.setBases(bases);
        read.setBaseQualities(quals);
        return read;
    }

    private void maskWindow(final byte[] baseCounts, final int windowStart, final int windowEnd, final byte[] bases, final byte[] quals) {
        if (baseCounts['A'] + baseCounts['T'] + baseCounts['N'] >= threshAT || baseCounts['G'] + baseCounts['C'] + baseCounts['N'] >= threshGC) {
            for (int j = windowStart; j <= windowEnd; j++) {
                bases[j] = 'N';
                quals[j] = MASKED_BASE_QUAL;
            }
        }
    }

}