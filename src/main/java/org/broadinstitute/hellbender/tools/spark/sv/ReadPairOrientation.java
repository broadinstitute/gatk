package org.broadinstitute.hellbender.tools.spark.sv;

import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVFastqUtils;
import org.broadinstitute.hellbender.tools.spark.sv.utils.Strand;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.util.Objects;

/**
 * Enumeration with all possible read-pair orientations.
 * <p>
 *     Notice the presence of the special item {@link #XX} to indicated an unknown or undetermined orientation
 *     (e.g. when one in the pair is unmapped).
 * </p>
 */
public enum ReadPairOrientation {

    /**
     * Left-Right orientation, that is, with Illumina reads both mates were sequenced pointing to each other which
     * is the expected orientation when the original templates co-aligns withe the reference.
     */
    LR,

    /**
     * Left-left orientation, that is, with Illumina<sup>&reg;</sup> pair-sequencing reads both mates were sequence pointing to the left.
     */
    LL,

    /**
     * Right-right orientation, that is, with Illumina<sup>&reg;</sup> pair-sequencing reads both mates were sequence pointing to the right.
     */
    RR,

    /**
     * Right-left orientation, that is, with Illumina<sup>&reg;</sup> pair-sequencing reads both mates were sequenced pointing away from each other.
     */
    RL,

    /**
     * Special enum element to represent an unknown read pair orientation.
     */
    XX;

    public static final ReadPairOrientation PROPER = LR;

    /**
     * Checks whether the read pair orientation indicates that the original template coaligns with the reference.
     * @return {@code true} iff so.
     */
    public boolean isProper() {
        return this == PROPER;
    }

    /**
     * Checks whether this read pair orientation indicates that the original template has an inversion with respect
     * to the reference.
     * @return {@code true} iff so.
     */
    public boolean impliesInversion() {
        return this == LL || this == RR;
    }

    /**
     * Checks whether this read pair orientation indicates that the original template contains either a (tandem) duplication or
     * or short range translocation with respect to the reference; it is impossible to distinguish between these two
     * cases with the pair orientation alone.
     * @return {@code true} iff so.
     */
    public boolean impliesTranslocationOrDuplication() {
        return this == RL;
    }

    /**
     * Checks whether this read pair orientation indicates a concrete defined orientation.
     * @return {@code true} iff so.
     */
    public boolean isDefined() {
        return this != XX;
    }

    /**
     * Checkes whether this read pair orientation instance indicates that the actual orientation is undefined or unknown.
     * @return {@code true} iff so.
     */
    public boolean isUndefined() {
        return this == XX;
    }

    /**
     * Figures out the pair orientation looking into a single read.
     * <p>
     *     This method relies on the mate mapping information present in the input read record
     *     to determine the relative positioning and orientation for the read's mate.
     * </p>
     * <p>
     *     If the paired flag is not set in the input read, this method returns the undefined
     *     orientation {@link #XX}. Needless to say, this is also the case iff at least one of the
     *     mates is flagged as unmapped or their are mapped to different reference sequences.
     * </p>
     *
     * @throws IllegalArgumentException if {@code r1} is {@code null}.
     * @param r1 the only input read.
     * @return never {@code null}.
     */
    public static ReadPairOrientation fromRead(final SAMRecord r1) {
        Utils.nonNull(r1);
        if (!r1.getReadPairedFlag() || r1.getReadUnmappedFlag() || r1.getMateUnmappedFlag()) {
            return XX;
        } else {
            if (!readAndMateMapToTheSameContig(r1)) {
                return XX;
            } else if (r1.getStart() <= r1.getMateAlignmentStart()) {
                return fromReadPairNegativeStrandFlags(r1.getReadNegativeStrandFlag(), r1.getMateNegativeStrandFlag());
            } else {
                return fromReadPairNegativeStrandFlags(r1.getMateNegativeStrandFlag(), r1.getReadNegativeStrandFlag());
            }
        }
    }

    /**
     * Returns the orientation of a coordinate sorted pair that map to the same contig given their strands.
     * <p>
     *     A {@code null} strand indicates that the strand information for that read in the pair is unknown.
     * </p>
     * @param left the member of the pair that maps left-most on that contig.
     * @param right the member of the pair that maps right-most on that contig.
     * @return never {@code null}, but {@link #XX} if the strand information is missing for any in the pair.
     */
    public static ReadPairOrientation fromStrands(final Strand left, final Strand right) {
        if (left == null || right == null) {
            return XX;
        } else {
            return fromReadPairNegativeStrandFlags(left == Strand.NEGATIVE,
                    right == Strand.NEGATIVE);
        }
    }

    private static boolean readAndMateMapToTheSameContig(final SAMRecord read) {
        final Integer readReferenceIndex = read.getReferenceIndex();
        final Integer mateReferenceIndex = read.getMateReferenceIndex();
        final boolean sameReference;
        if (readReferenceIndex == null || mateReferenceIndex == null)  {
            sameReference = Objects.equals(read.getReferenceName(), read.getMateReferenceName());
        } else {
            sameReference = readReferenceIndex.equals(mateReferenceIndex);
        }
        return sameReference;
    }

    private static ReadPairOrientation fromReadPairNegativeStrandFlags(final boolean r1NegativeStrandFlag,
                                                                       final boolean r2NegativeStrandFlag) {
        if (!r1NegativeStrandFlag) {
            return r2NegativeStrandFlag ? LR : LL;
        } else {
            return r2NegativeStrandFlag ? RR : RL;
        }
    }
}
