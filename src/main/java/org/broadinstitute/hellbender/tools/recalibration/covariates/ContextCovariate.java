package org.broadinstitute.hellbender.tools.recalibration.covariates;

import htsjdk.samtools.SAMRecord;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.recalibration.ReadCovariates;
import org.broadinstitute.hellbender.tools.recalibration.RecalibrationArgumentCollection;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.clipping.ClippingRepresentation;
import org.broadinstitute.hellbender.utils.clipping.ReadClipper;

import java.util.ArrayList;
import java.util.Arrays;

public final class ContextCovariate implements Covariate {
    private final static Logger logger = LogManager.getLogger(ContextCovariate.class);

    private int mismatchesContextSize;
    private int indelsContextSize;

    private int mismatchesKeyMask;
    private int indelsKeyMask;

    private static final int LENGTH_BITS = 4;
    private static final int LENGTH_MASK = 15;

    // the maximum context size (number of bases) permitted; we need to keep the leftmost base free so that values are
    // not negative and we reserve 4 more bits to represent the length of the context; it takes 2 bits to encode one base.
    static final private int MAX_DNA_CONTEXT = 13;
    private byte LOW_QUAL_TAIL;

    // Initialize any member variables using the command-line arguments passed to the walkers
    @Override
    public void initialize(final RecalibrationArgumentCollection RAC) {
        mismatchesContextSize = RAC.MISMATCHES_CONTEXT_SIZE;
        indelsContextSize = RAC.INDELS_CONTEXT_SIZE;

        logger.info("\t\tContext sizes: base substitution model " + mismatchesContextSize + ", indel substitution model " + indelsContextSize);

        if (mismatchesContextSize > MAX_DNA_CONTEXT)
            throw new UserException.BadArgumentValue("mismatches_context_size", String.format("context size cannot be bigger than %d, but was %d", MAX_DNA_CONTEXT, mismatchesContextSize));
        if (indelsContextSize > MAX_DNA_CONTEXT)
            throw new UserException.BadArgumentValue("indels_context_size", String.format("context size cannot be bigger than %d, but was %d", MAX_DNA_CONTEXT, indelsContextSize));

        LOW_QUAL_TAIL = RAC.LOW_QUAL_TAIL;
        
        if (mismatchesContextSize <= 0 || indelsContextSize <= 0)
            throw new UserException(String.format("Context size must be positive, if you don't want to use the context covariate, just turn it off instead. Mismatches: %d Indels: %d", mismatchesContextSize, indelsContextSize));

        mismatchesKeyMask = createMask(mismatchesContextSize);
        indelsKeyMask = createMask(indelsContextSize);
    }

    @Override
    public void recordValues(final SAMRecord read, final ReadCovariates values) {

        // store the original bases and then write Ns over low quality ones
        final byte[] originalBases = Arrays.copyOf(read.getReadBases(), read.getReadBases().length);
        // Write N's over the low quality tail of the reads to avoid adding them into the context
        final SAMRecord clippedRead = ReadClipper.clipLowQualEnds(read, LOW_QUAL_TAIL, ClippingRepresentation.WRITE_NS);
        
        final boolean negativeStrand = clippedRead.getReadNegativeStrandFlag();
        byte[] bases = clippedRead.getReadBases();
        if (negativeStrand)
            bases = BaseUtils.simpleReverseComplement(bases);

        final ArrayList<Integer> mismatchKeys = contextWith(bases, mismatchesContextSize, mismatchesKeyMask);
        final ArrayList<Integer> indelKeys = contextWith(bases, indelsContextSize, indelsKeyMask);

        final int readLength = bases.length;

        // this is necessary to ensure that we don't keep historical data in the ReadCovariates values
        // since the context covariate may not span the entire set of values in read covariates
        // due to the clipping of the low quality bases
        if ( readLength != originalBases.length ) {
            // don't both zeroing out if we are going to overwrite the whole array
            for ( int i = 0; i < originalBases.length; i++ )
                // this base has been clipped off, so zero out the covariate values here
                values.addCovariate(0, 0, 0, i);
        }

        for (int i = 0; i < readLength; i++) {
            final int readOffset = (negativeStrand ? readLength - i - 1 : i);
            final int indelKey = indelKeys.get(i);
            values.addCovariate(mismatchKeys.get(i), indelKey, indelKey, readOffset);
        }

        // put the original bases back in
        read.setReadBases(originalBases);
    }

    // Used to get the covariate's value from input csv file during on-the-fly recalibration
    @Override
    public final Object getValue(final String str) {
        return str;
    }

    @Override
    public String formatKey(final int key) {
        if (key == -1) // this can only happen in test routines because we do not propagate null keys to the csv file
            return null;

        return contextFromKey(key);
    }

    @Override
    public int keyFromValue(final Object value) {
        return keyFromContext((String) value);
    }

    private static int createMask(final int contextSize) {
        int mask = 0;
        // create 2*contextSize worth of bits
        for (int i = 0; i < contextSize; i++)
            mask = (mask << 2) | 3;
        // shift 4 bits to mask out the bits used to encode the length
        return mask << LENGTH_BITS;
    }

    /**
     * calculates the context of a base independent of the covariate mode (mismatch, insertion or deletion)
     *
     * @param bases       the bases in the read to build the context from
     * @param contextSize context size to use building the context
     * @param mask        mask for pulling out just the context bits
     */
    private static ArrayList<Integer> contextWith(final byte[] bases, final int contextSize, final int mask) {

        final int readLength = bases.length;
        final ArrayList<Integer> keys = new ArrayList<Integer>(readLength);

        // the first contextSize-1 bases will not have enough previous context
        for (int i = 1; i < contextSize && i <= readLength; i++)
            keys.add(-1);

        if (readLength < contextSize)
            return keys;

        final int newBaseOffset = 2 * (contextSize - 1) + LENGTH_BITS;

        // get (and add) the key for the context starting at the first base
        int currentKey = keyFromContext(bases, 0, contextSize);
        keys.add(currentKey);

        // if the first key was -1 then there was an N in the context; figure out how many more consecutive contexts it affects
        int currentNPenalty = 0;
        if (currentKey == -1) {
            currentKey = 0;
            currentNPenalty = contextSize - 1;
            int offset = newBaseOffset;
            while (bases[currentNPenalty] != 'N') {
                final int baseIndex = BaseUtils.simpleBaseToBaseIndex(bases[currentNPenalty]);
                currentKey |= (baseIndex << offset);
                offset -= 2;
                currentNPenalty--;
            }
        }

        for (int currentIndex = contextSize; currentIndex < readLength; currentIndex++) {
            final int baseIndex = BaseUtils.simpleBaseToBaseIndex(bases[currentIndex]);
            if (baseIndex == -1) { // ignore non-ACGT bases
                currentNPenalty = contextSize;
                currentKey = 0; // reset the key
            } else {
                // push this base's contribution onto the key: shift everything 2 bits, mask out the non-context bits, and add the new base and the length in
                currentKey = (currentKey >> 2) & mask;
                currentKey |= (baseIndex << newBaseOffset);
                currentKey |= contextSize;
            }

            if (currentNPenalty == 0) {
                keys.add(currentKey);
            } else {
                currentNPenalty--;
                keys.add(-1);
            }
        }

        return keys;
    }

    public static int keyFromContext(final String dna) {
        return keyFromContext(dna.getBytes(), 0, dna.length());
    }

    /**
     * Creates a int representation of a given dna string.
     *
     * @param dna    the dna sequence
     * @param start  the start position in the byte array (inclusive)
     * @param end    the end position in the array (exclusive)
     * @return the key representing the dna sequence
     */
    private static int keyFromContext(final byte[] dna, final int start, final int end) {

        int key = end - start;
        int bitOffset = LENGTH_BITS;
        for (int i = start; i < end; i++) {
            final int baseIndex = BaseUtils.simpleBaseToBaseIndex(dna[i]);
            if (baseIndex == -1) // ignore non-ACGT bases
                return -1;
            key |= (baseIndex << bitOffset);
            bitOffset += 2;
        }
        return key;
    }

    /**
     * Converts a key into the dna string representation.
     *
     * @param key    the key representing the dna sequence
     * @return the dna sequence represented by the key
     */
    public static String contextFromKey(final int key) {
        if (key < 0)
            throw new GATKException("dna conversion cannot handle negative numbers. Possible overflow?");

        final int length = key & LENGTH_MASK; // the first bits represent the length (in bp) of the context
        int mask = 48; // use the mask to pull out bases
        int offset = LENGTH_BITS;

        StringBuilder dna = new StringBuilder();
        for (int i = 0; i < length; i++) {
            final int baseIndex = (key & mask) >> offset;
            dna.append((char)BaseUtils.baseIndexToSimpleBase(baseIndex));
            mask = mask << 2; // move the mask over to the next 2 bits
            offset += 2;
        }

        return dna.toString();
    }

    @Override
    public int maximumKeyValue() {
        // the maximum value is T (11 in binary) for each base in the context
        int length = Math.max(mismatchesContextSize, indelsContextSize);  // the length of the context
        int key = length;
        int bitOffset = LENGTH_BITS;
        for (int i = 0; i <length ; i++) {
            key |= (3 << bitOffset);
            bitOffset += 2;
        }
        return key;
    }
}
