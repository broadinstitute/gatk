package org.broadinstitute.hellbender.utils.recalibration.covariates;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntList;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.clipping.ClippingRepresentation;
import org.broadinstitute.hellbender.utils.clipping.ReadClipper;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.recalibration.RecalibrationArgumentCollection;

import java.util.Arrays;

public class ContextCovariate implements Covariate {
    private static final long serialVersionUID = 1L;
    private static final Logger logger = LogManager.getLogger(ContextCovariate.class);

    private final int mismatchesContextSize;
    private final int indelsContextSize;

    private final int mismatchesKeyMask;
    private final int indelsKeyMask;

    protected static final int LENGTH_BITS = 4;
    private static final int LENGTH_MASK = 15;


    // the maximum context size (number of bases) permitted; we need to keep the leftmost base free so that values are
    // not negative and we reserve 4 more bits to represent the length of the context; it takes 2 bits to encode one base.
    private static final int MAX_DNA_CONTEXT = 13;
    private final byte lowQualTail;

    /**
     * amount of lookhead allocated from the context size. The first by of the lookahead is always taken from the reference
     *
     * For example, if the context size is 3 and the amount of lookahead allocated is 2 then
     * the context will be composed of 1 byte before the location, 1 byte from the reference and 1 byte after the location
     */
    private int lookaheadSize;
    private ReferenceDataSource referenceDataSource;
    private boolean altEmbedded;

    public ContextCovariate(final RecalibrationArgumentCollection RAC){
        mismatchesContextSize = RAC.MISMATCHES_CONTEXT_SIZE;
        indelsContextSize = RAC.INDELS_CONTEXT_SIZE;
        logger.debug("\t\tContext sizes: base substitution model " + mismatchesContextSize + ", indel substitution model " + indelsContextSize);

        if (mismatchesContextSize > MAX_DNA_CONTEXT) {
            throw new CommandLineException.BadArgumentValue("mismatches_context_size", String.format("context size cannot be bigger than %d, but was %d", MAX_DNA_CONTEXT, mismatchesContextSize));
        }
        if (indelsContextSize > MAX_DNA_CONTEXT) {
            throw new CommandLineException.BadArgumentValue("indels_context_size", String.format("context size cannot be bigger than %d, but was %d", MAX_DNA_CONTEXT, indelsContextSize));
        }

        lowQualTail = RAC.LOW_QUAL_TAIL;

        if (mismatchesContextSize <= 0 || indelsContextSize <= 0) {
            throw new CommandLineException(String.format("Context size must be positive. Mismatches: %d Indels: %d", mismatchesContextSize, indelsContextSize));
        }

        mismatchesKeyMask = createMask(mismatchesContextSize);
        indelsKeyMask = createMask(indelsContextSize);

        if ( RAC.EXTENDED_CONTEXT_LOOKAHEAD != 0 ) {
            Utils.nonNull(RAC.EXTENDED_CONTEXT_REFERENCE, "extended context reference can not be null");
            Utils.validate(RAC.EXTENDED_CONTEXT_LOOKAHEAD > 0, "lookahead must be positive");
            Utils.validate(RAC.EXTENDED_CONTEXT_LOOKAHEAD <= Math.min(RAC.MISMATCHES_CONTEXT_SIZE, RAC.INDELS_CONTEXT_SIZE), "lookahead can not be larger than the context");

            this.lookaheadSize = RAC.EXTENDED_CONTEXT_LOOKAHEAD;
            this.referenceDataSource = ReferenceDataSource.of(RAC.EXTENDED_CONTEXT_REFERENCE.toPath());
            this.altEmbedded = RAC.EXTENDED_CONTEXT_ALT_EMBEDDED;
        }
    }

    @Override
    public void recordValues(final GATKRead read, final SAMFileHeader header, final ReadCovariates values, final boolean recordIndelValues) {

        final int originalReadLength = read.getLength();

        // store the original bases and then write Ns over low quality ones
        final byte[] strandedClippedBases = getStrandedClippedBytes(read, lowQualTail);  //Note: this makes a copy of the read
        final byte[] refBases = isExtended() ? getReadReferenceBases(read) : null;

        //Note: we're using a non-standard library here because boxing came up on profiling as taking 20% of time in applyBQSR.
        //IntList avoids boxing
        final IntList mismatchKeys = isExtended()
                ? contextWith(strandedClippedBases, mismatchesContextSize, mismatchesKeyMask, refBases)
                : contextWith(strandedClippedBases, mismatchesContextSize, mismatchesKeyMask);

        final int readLengthAfterClipping = strandedClippedBases.length;

        // this is necessary to ensure that we don't keep historical data in the ReadCovariates values
        // since the context covariate may not span the entire set of values in read covariates
        // due to the clipping of the low quality bases
        if ( readLengthAfterClipping != originalReadLength) {
            // don't bother zeroing out if we are going to overwrite the whole array
            for ( int i = 0; i < originalReadLength; i++ ){
                // this base has been clipped off, so zero out the covariate values here
                values.addCovariate(0, 0, 0, i);
            }
        }

        final boolean negativeStrand = read.isReverseStrand();

        //Note: duplicated the loop to avoid checking recordIndelValues on each iteration
        if (recordIndelValues) {
            final IntList indelKeys = isExtended()
                    ? contextWith(strandedClippedBases, indelsContextSize, indelsKeyMask, refBases)
                    : contextWith(strandedClippedBases, indelsContextSize, indelsKeyMask);
            for (int i = 0; i < readLengthAfterClipping; i++) {
                final int readOffset = getStrandedOffset(negativeStrand, i, readLengthAfterClipping);
                final int indelKey = indelKeys.getInt(i);
                values.addCovariate(mismatchKeys.getInt(i), indelKey, indelKey, readOffset);
            }
        } else {
            for (int i = 0; i < readLengthAfterClipping; i++) {
                final int readOffset = getStrandedOffset(negativeStrand, i, readLengthAfterClipping);
                values.addCovariate(mismatchKeys.getInt(i), 0, 0, readOffset);
            }
        }
    }

    /**
     * Helper method: computes the correct offset to use in computations of covariate values.
     * @param isNegativeStrand is the read on the negative strand
     * @param offset 0-based index of the base in the read
     * @param readLength length of the read
     * @return
     */
    public static int getStrandedOffset(final boolean isNegativeStrand, final int offset, final int readLength) {
        return isNegativeStrand ? (readLength - offset - 1) : offset;
    }

    /**
     * Given a read, clips low quality ends (by overwriting with N) and returns the underlying bases, after
     * reverse-complementing for negative-strand reads.
     * @param read the read
     * @param lowQTail every base quality lower than or equal to this in the tail of the read will be replaced with N.
     * @return bases of the read (Could be an empty array if all bases are below lowQTail).
     */
    @VisibleForTesting
    static byte[] getStrandedClippedBytes(final GATKRead read, final byte lowQTail) {

        // Write N's over the low quality tail of the reads to avoid adding them into the context
        final GATKRead clippedRead = ReadClipper.clipLowQualEnds(read, lowQTail, ClippingRepresentation.WRITE_NS);

        final byte[] bases = clippedRead.getBases();
        if (read.isReverseStrand()) {
            return BaseUtils.simpleReverseComplement(bases);
        } else {
            return bases;
        }
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
        for (int i = 0; i < contextSize; i++) {
            mask = (mask << 2) | 3;
        }
        // shift 4 bits to mask out the bits used to encode the length
        return mask << LENGTH_BITS;
    }

    /**
     * calculates the context of a base independent of the covariate mode (mismatch, insertion or deletion)
     *  @param bases       the bases in the read to build the context from
     * @param contextSize context size to use building the context
     * @param mask        mask for pulling out just the context bits
     */
    protected IntList contextWith(final byte[] bases, final int contextSize, final int mask) {

        final int readLength = bases.length;

        //Note: we use a specialized collection to avoid the cost of boxing and unboxing that otherwise comes up on the profiler.
        final IntList keys = new IntArrayList(readLength);

        // the first contextSize-1 bases will not have enough previous context
        for (int i = 1; i < contextSize && i <= readLength; i++) {
            keys.add(-1);
        }

        if (readLength < contextSize) {
            return keys;
        }

        final int newBaseOffset = 2 * (contextSize - 1) + LENGTH_BITS;

        // get (and add) the key for the context starting at the first base
        int currentKey = keyFromContext(bases, 0, contextSize);
        keys.add(currentKey);

        // if the first key was -1 then there was an non-ACGT in the context; figure out how many more consecutive contexts it affects
        int currentNPenalty = 0;
        if (currentKey == -1) {
            currentKey = 0;
            currentNPenalty = contextSize - 1;
            int offset = newBaseOffset;
            int baseIndex;
            while ((baseIndex = BaseUtils.simpleBaseToBaseIndex(bases[currentNPenalty])) != -1) {
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

    public int keyFromContext(final String dna) {
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
            if (baseIndex == -1) { // ignore non-ACGT bases
                return -1;
            }
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
        if (key < 0) {
            throw new GATKException("dna conversion cannot handle negative numbers. Possible overflow?");
        }

        final int length = key & LENGTH_MASK; // the first bits represent the length (in bp) of the context
        int mask = 48; // use the mask to pull out bases
        int offset = LENGTH_BITS;

        final StringBuilder dna = new StringBuilder(length);
        for (int i = 0; i < length; i++) {
            final int baseIndex = (key & mask) >> offset;
            dna.append((char)BaseUtils.baseIndexToSimpleBase(baseIndex));
            mask <<= 2; // move the mask over to the next 2 bits
            offset += 2;
        }

        return dna.toString();
    }

    @Override
    public int maximumKeyValue() {
        // the maximum value is T (11 in binary) for each base in the context
        final int length = Math.max(mismatchesContextSize, indelsContextSize)
                + (altEmbedded ? 1 : 0);  // the length of the context
        int key = length;
        int bitOffset = LENGTH_BITS;
        for (int i = 0; i <length ; i++) {
            key |= (3 << bitOffset);
            bitOffset += 2;
        }
        return key;
    }

    private boolean isExtended() {
        return lookaheadSize != 0;
    }

    /**
     * calculates the context of a base independent of the covariate mode (mismatch, insertion or deletion)
     *  @param bases       the bases in the read to build the context from
     * @param contextSize context size to use building the context
     * @param mask        mask for pulling out just the context bits
     */
    protected IntList contextWith(final byte[] bases, final int contextSize, final int mask, final byte[] refBases) {

        if ( bases.length == 0 ) {
            return new IntArrayList(0);
        }

        Utils.validate(bases.length == refBases.length, "read and references bases array not same length");

        final int readLength = bases.length;

        //Note: we use a specialized collection to avoid the cost of boxing and unboxing that otherwise comes up on the profiler.
        final IntList keys = new IntArrayList(readLength);

        // determine where context will actually start and end in relationship to the base in question
        if ( readLength < contextSize ) {
            for ( int i = 0 ; i < readLength ; i++ ) {
                keys.add(-1);
            }
            return keys;
        }
        final int contextBefore =  contextSize - lookaheadSize;   // how many bases before current base are part of the context
        final int contextAfter = lookaheadSize - 1;              // how many bases after current base are part of the context

        // before context can be established
        for ( int i = 0 ; i < contextBefore ; i++ ) {
            keys.add(-1);
        }

        // main loop, for now compute for each position
        for ( int i = contextBefore ; i < readLength - contextAfter ; i++ ) {

            int currentKey = keyFromContext(bases, i - contextBefore, i + contextAfter + 1, i, refBases[i]);
            keys.add(currentKey);
        }

        // after context can be established
        for ( int i = readLength - contextAfter ; i < readLength ; i++ ) {
            keys.add(-1);
        }

        Utils.validate(keys.size() == readLength, "generated keys and readLength differ");;

        return keys;
    }

    protected byte[] getReadReferenceBases(GATKRead read) {
        final byte[]    refBases = referenceDataSource.queryAndPrefetch(read.getContig(), read.getStart(), read.getEnd()).getBases();
        final byte[]    readRefBases = new byte[read.getLength()];
        int             readOfs = 0;
        int             refOfs = 0;
        for (CigarElement elem : read.getCigarElements() ) {
            final CigarOperator op = elem.getOperator();
            final int length = elem.getLength();
            if ( op.consumesReadBases() && op.consumesReferenceBases() ) {
                // simple case - copy from reference
                System.arraycopy(refBases, refOfs, readRefBases, readOfs, length);
                refOfs += length;
                readOfs += length;
            } else if ( op.consumesReadBases() ) {
                // has read bases but not reference bases -> INS. supplement reference with N
                Arrays.fill(readRefBases, readOfs, readOfs + length, (byte)'N');
                readOfs += length;
            } else if ( op.consumesReferenceBases() ) {
                // has reference but not read bases -> DEL. skip on reference
                refOfs += length;
            }
        }

        Utils.validate(readOfs == read.getLength(), "did not read end of read");
        Utils.validate(refOfs == refBases.length, "did not reach end of reference");

        // reverse complement?
        if ( read.isReverseStrand() ) {
            return BaseUtils.simpleReverseComplement(readRefBases);
        } else {
            return readRefBases;
        }
    }

    private int keyFromContext(final byte[] dna, final int start, final int end, final int refIndex, final byte refBase) {

        int key = end - start + (altEmbedded ? 1 : 0);
        int bitOffset = LENGTH_BITS;
        for (int i = start; i < end; i++) {
            final int baseIndex = BaseUtils.simpleBaseToBaseIndex((i != refIndex) ? dna[i] : refBase);
            if (baseIndex == -1) { // ignore non-ACGT bases
                return -1;
            }
            key |= (baseIndex << bitOffset);
            bitOffset += 2;

            if ( altEmbedded && i == refIndex ) {
                final int altBaseIndex = BaseUtils.simpleBaseToBaseIndex(dna[i]);
                if (altBaseIndex == -1) { // ignore non-ACGT bases
                    return -1;
                }
                key |= (altBaseIndex << bitOffset);
                bitOffset += 2;
            }
        }
        return key;
    }
}
