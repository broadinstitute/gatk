package org.broadinstitute.hellbender.transformers;

import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.clipping.ClippingOp;
import org.broadinstitute.hellbender.utils.clipping.ClippingRepresentation;
import org.broadinstitute.hellbender.utils.clipping.ReadClipper;
import org.broadinstitute.hellbender.utils.read.GATKRead;

/**
 * Masks read bases with a supra-threshold number of A/T's or G/C's within a given window size.
 */
public final class SimpleRepeatMaskTransformer implements ReadTransformer {
    private static final long serialVersionUID = 1L;

    private final int windowSize;
    private final int threshAT;
    private final int threshGC;

    public SimpleRepeatMaskTransformer(final int threshAT, final int threshGC, final int windowSize) {
        Utils.validateArg(threshAT <= windowSize, "SimpleRepeatMaskTransformer threshold AT greater than window size.");
        Utils.validateArg(threshGC <= windowSize, "SimpleRepeatMaskTransformer threshold GC greater than window size.");
        this.threshAT = threshAT;
        this.threshGC = threshGC;
        this.windowSize = windowSize;
    }

    @Override
    public GATKRead apply(GATKRead read) {
        final byte[] originalBases = read.getBases().clone();
        byte[] baseCounts = new byte[255];

        //Initialize with the first window
        int windowStart = 0;
        int windowEnd = Math.min(windowSize, originalBases.length);
        for (int i = windowStart; i < windowEnd; i++) {
            baseCounts[originalBases[i]]++;
        }
        int maskStart = 0;
        int maskEnd = 0;
        if (isSimpleRepeat(baseCounts)) {
            maskEnd = windowEnd;
        }

        //March the window through the read, one base at a time
        while (windowEnd < originalBases.length) {
            baseCounts[originalBases[windowStart++]]--; //Subtract old first base
            baseCounts[originalBases[windowEnd++]]++; //Add new final base
            if (isSimpleRepeat(baseCounts)) {
                if (windowEnd == maskEnd + 1) {
                    maskEnd++;
                } else {
                    read = maskRead(read, maskStart, maskEnd);
                    maskStart = windowStart;
                    maskEnd = windowEnd;
                }
            }
        }
        read = maskRead(read, maskStart, maskEnd);
        return read;
    }

    private boolean isSimpleRepeat(final byte[] baseCounts) {
        final int countN = baseCounts['N'] + baseCounts['n'];
        final int countATN = baseCounts['A'] + baseCounts['a'] + baseCounts['T'] + baseCounts['t'] + countN;
        final int countGCN = baseCounts['G'] + baseCounts['g'] + baseCounts['C'] + baseCounts['c'] + countN;
        return countATN >= threshAT || countGCN >= threshGC;
    }

    private GATKRead maskRead(final GATKRead read, final int maskStart, final int maskEnd) {
        final ReadClipper readClipper = new ReadClipper(read);
        readClipper.addOp(new ClippingOp(maskStart, maskEnd - 1));
        return readClipper.clipRead(ClippingRepresentation.WRITE_NS_Q0S);
    }

}