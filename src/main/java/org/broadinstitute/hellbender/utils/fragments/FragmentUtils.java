package org.broadinstitute.hellbender.utils.fragments;

import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.util.QualityUtil;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.util.OptionalInt;

public final class FragmentUtils {
    private FragmentUtils() {}

    public final static double DEFAULT_PCR_SNV_ERROR_RATE = 1e-4;
    public final static int DEFAULT_PCR_SNV_ERROR_QUAL = QualityUtil.getPhredScoreFromErrorProbability(DEFAULT_PCR_SNV_ERROR_RATE);
    public final static int HALF_OF_DEFAULT_PCR_SNV_ERROR_QUAL = DEFAULT_PCR_SNV_ERROR_QUAL / 2;

    /**
     * Fix two overlapping reads from the same fragment by adjusting base qualities, if possible
     *
     *  Looks at the bases and alignment, and tries its best to create adjusted base qualities so that the observations
     * are not treated independently.  Sets the qualities of firstRead and secondRead to mimic a merged read or
     * nothing if the algorithm cannot create a meaningful one
     * @param pair two overlapping paired reads
     * @param setConflictingToZero if true, set base qualities to zero when mates have different base at overlapping position
     * @param halfOfPcrSnvQual half of phred-scaled quality of substitution errors from PCR. May not be negative.
     * @param halfOfPcrIndelQual half of phred-scaled quality of indel errors from PCR. May not be negative.
     */
    public static void adjustQualsOfOverlappingPairedFragments(final Pair<GATKRead, GATKRead> pair, final boolean setConflictingToZero,
                                                               final OptionalInt halfOfPcrSnvQual, final OptionalInt   halfOfPcrIndelQual) {
        final boolean inOrder = pair.getLeft().getSoftStart() < pair.getRight().getSoftStart();
        final GATKRead firstRead = inOrder ? pair.getLeft() : pair.getRight();
        final GATKRead secondRead = inOrder ? pair.getRight() : pair.getLeft();

        Utils.nonNull(firstRead);
        Utils.nonNull(secondRead);
        Utils.validateArg(firstRead.getName().equals(secondRead.getName()), () ->
                "attempting to merge two reads with different names " + firstRead + " and " + secondRead);

        // don't adjust fragments that do not overlap
        if (firstRead.getEnd() < secondRead.getStart() || !firstRead.getContig().equals(secondRead.getContig())) {
            return;
        }

        // the offset and cigar operator in the first read at the start of the left read
        final Pair<Integer, CigarOperator> offsetAndOperator = ReadUtils.getReadIndexForReferenceCoordinate(firstRead, secondRead.getStart());
        final CigarOperator operator = offsetAndOperator.getRight();
        final int offset = offsetAndOperator.getLeft();
        if (offset == ReadUtils.READ_INDEX_NOT_FOUND || operator.isClipping()) { // no overlap or only overlap in clipped region
            return;
        }

        // Compute the final aligned base indexes for both since there might be right base softclips
        final int firstReadEndBase = ReadUtils.getReadIndexForReferenceCoordinate(firstRead, firstRead.getEnd()).getLeft();
        final int secondReadEndBase = ReadUtils.getReadIndexForReferenceCoordinate(secondRead, secondRead.getEnd()).getLeft();

        // TODO: we should be careful about the case where {@code operator} is a deletion; that is, when the second read start falls in a deletion of the first read
        // TODO: however, the issue is bigger than simply getting the start correctly, because the code below assumes that all bases of both reads are aligned in their overlap.
        // TODO: Any indel that occurs in one read and not the other will spoil things.  Really, the correct thing to do is a Smith-Waterman (or other) alignment of the reads
        // TODO: in their overlap and correct the double-counting for all aligned bases.
        // TODO: a cheaper solution would be to cap all quals in the overlap region by half of the PCR qual.
        final int firstReadStop = offset;
        final int secondOffset = ReadUtils.getReadIndexForReferenceCoordinate(secondRead, secondRead.getStart()).getLeft(); //This operation handles softclipped bases in the qual/base array
        final int numOverlappingBases = Math.min(firstReadEndBase - firstReadStop, secondReadEndBase - secondOffset) + 1; // Add 1 here because if R1 ends on the same base that R2 starts then there is 1 base of overlap not 0

        final byte[] firstReadBases = firstRead.getBases();
        final byte[] firstReadQuals = firstRead.getBaseQualities();
        final byte[] secondReadBases = secondRead.getBases();
        final byte[] secondReadQuals = secondRead.getBaseQualities();

        final int halfOfPcrErrorQual = halfOfPcrSnvQual.orElse(HALF_OF_DEFAULT_PCR_SNV_ERROR_QUAL);

        for (int i = 0; i < numOverlappingBases; i++) {

            final int firstReadIndex = firstReadStop + i;
            final int secondReadIndex = secondOffset + i;
            final byte firstReadBase = firstReadBases[firstReadIndex];
            final byte secondReadBase = secondReadBases[secondReadIndex];

            if (firstReadBase == secondReadBase) {
                firstReadQuals[firstReadIndex] = (byte) Math.min(firstReadQuals[firstReadIndex], halfOfPcrErrorQual);
                secondReadQuals[secondReadIndex] = (byte) Math.min(secondReadQuals[secondReadIndex], halfOfPcrErrorQual);
            } else if (setConflictingToZero) {
                // If downstream processing forces read pairs to support the same haplotype, setConflictingToZero should be false
                // because the original base qualities of conflicting bases, when pegged to the same haplotype, will
                // automatically weaken the strength of one another's evidence.  Furthermore, if one base if low quality
                // and one is high it will essentially ignore the low quality base without compromising the high-quality base
                firstReadQuals[firstReadIndex] = 0;
                secondReadQuals[secondReadIndex] = 0;
            }
        }
        firstRead.setBaseQualities(firstReadQuals);
        secondRead.setBaseQualities(secondReadQuals);

        if (halfOfPcrIndelQual.isPresent()) {
            final int maxIndelQual = halfOfPcrIndelQual.getAsInt();
            final byte[] firstReadInsertionQuals = ReadUtils.getBaseInsertionQualities(firstRead);
            final byte[] firstReadDeletionQuals = ReadUtils.getBaseDeletionQualities(firstRead);
            final byte[] secondReadInsertionQuals = ReadUtils.getBaseInsertionQualities(secondRead);
            final byte[] secondReadDeletionQuals = ReadUtils.getBaseDeletionQualities(secondRead);

            for (int i = 0; i < numOverlappingBases; i++) {
                final int firstReadIndex = firstReadStop + i;
                final int secondReadIndex = secondOffset + i;
                firstReadDeletionQuals[firstReadIndex] = (byte) Math.min(firstReadDeletionQuals[firstReadIndex], maxIndelQual);
                firstReadInsertionQuals[firstReadIndex] = (byte) Math.min(firstReadInsertionQuals[firstReadIndex], maxIndelQual);
                secondReadDeletionQuals[secondReadIndex] = (byte) Math.min(secondReadDeletionQuals[secondReadIndex], maxIndelQual);
                secondReadInsertionQuals[secondReadIndex] = (byte) Math.min(secondReadInsertionQuals[secondReadIndex], maxIndelQual);
            }

            ReadUtils.setDeletionBaseQualities(firstRead, firstReadDeletionQuals);
            ReadUtils.setInsertionBaseQualities(firstRead, firstReadInsertionQuals);
            ReadUtils.setDeletionBaseQualities(secondRead, secondReadDeletionQuals);
            ReadUtils.setInsertionBaseQualities(secondRead, secondReadInsertionQuals);
        }
    }
}
