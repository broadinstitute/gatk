package org.broadinstitute.hellbender.utils.pileup;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.fragments.FragmentCollection;

import java.util.*;
import java.util.function.Predicate;
import java.util.stream.Stream;

/**
 * Abstract class for tracking the pileup elements from different sources.
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public abstract class PileupElementTracker implements Iterable<PileupElement> {

    /** Constant used by samtools to downgrade a quality for overlapping reads that disagrees in their base. */
    public static final double SAMTOOLS_OVERLAP_LOW_CONFIDENCE = 0.8;

    /** Comparator for sorting by start read position*/
    protected static final Comparator<PileupElement> READ_START_COMPARATOR = (l, r) -> Integer.compare(l.getRead().getStart(), r.getRead().getStart());

    /**
     * Gets a stream over the elements in this tracker without any specific order.
     */
    public abstract Stream<PileupElement> getElementStream();

    /** Returns the number of elements in the tracker. */
    public final int size() {
        return (int) getElementStream().count();
    }

    /**
     * Iterates over the PileupElements without assuming any specific ordering.
     *
     * @return the iterator over the elements in this tracker.
     */
    @Override
    public final Iterator<PileupElement> iterator() {
        return getElementStream().iterator();
    }

    /**
     * Iterator over sorted by read start PileupElements.
     */
    public abstract Iterator<PileupElement> sortedIterator();

    /**
     * Gets a set of the samples represented in this pileup.
     * Note: contains null if a read has a null read group or a null sample name.
     */
    public abstract Set<String> getSamples(final SAMFileHeader header);

    /**
     * Splits the tracker by sample.
     *
     * @param header            the header to retrieve the samples from
     * @return an element tracker splitted by sample (including empty trackers for a sample)
     */
    public abstract PileupElementTracker splitBySample(final SAMFileHeader header);

    /**
     * Make a new tracker consisting of elements that satisfy the predicate.
     * NOTE: the new tracker will not be independent of the old one (no deep copy of the underlying data is performed).
     */
    public abstract PileupElementTracker makeFilteredTracker(final Predicate<PileupElement> filter);

    /**
     * Make a new tracker from elements whose reads belong to the given sample.
     *
     * Passing null sample as an argument retrieves reads whose read group or sample name is {@code null}
     * NOTE: the new tracker will not be independent of the old one (no deep copy of the underlying data is performed).
     */
    public abstract PileupElementTracker getTrackerForSample(final String sample, final SAMFileHeader header);

    /**
     * Fixes the quality of all the elements that come from an overlapping pair in the same way as
     * samtools does {@see tweak_overlap_quality function in
     * <a href="https://github.com/samtools/htslib/blob/master/sam.c">samtools</a>}.
     * <p>
     * Setting the quality of one of the bases to 0 effectively removes the redundant base for
     * calling. In addition, if the bases overlap we have increased confidence if they agree (or
     * reduced if they don't). Thus, the algorithm proceeds as following:
     * <p>
     * 1. If the bases are the same, the quality of the first element is the sum of both qualities
     * and the quality of the second is reduced to 0.
     * 2. If the bases are different, the base with the highest quality is reduced with a factor of
     * 0.8, and the quality of the lowest is reduced to 0.
     * <p>
     * Note: Resulting qualities higher than {@link QualityUtils#MAX_SAM_QUAL_SCORE} are capped.
     */
    public void fixOverlaps() {
        final FragmentCollection<PileupElement> fragments = FragmentCollection.create(this);
        fragments.getOverlappingPairs().stream()
                .forEach(
                        elements -> fixPairOverlappingQualities(elements.get(0), elements.get(1))
                );
    }

    /**
     * Fixes the quality of two elements that come from an overlapping pair in the same way as
     * samtools does {@see tweak_overlap_quality function in
     * <a href="https://github.com/samtools/htslib/blob/master/sam.c">samtools</a>}.
     * The only difference with the samtools API is the cap for high values ({@link QualityUtils#MAX_SAM_QUAL_SCORE}).
     */
    @VisibleForTesting
    static void fixPairOverlappingQualities(final PileupElement firstElement,
            final PileupElement secondElement) {
        // only if they do not represent deletions
        if (!secondElement.isDeletion() && !firstElement.isDeletion()) {
            final byte[] firstQuals = firstElement.getRead().getBaseQualities();
            final byte[] secondQuals = secondElement.getRead().getBaseQualities();
            if (firstElement.getBase() == secondElement.getBase()) {
                // if both have the same base, extra confidence in the firts of them
                firstQuals[firstElement.getOffset()] =
                        (byte) (firstQuals[firstElement.getOffset()] + secondQuals[secondElement
                                .getOffset()]);
                // cap to maximum byte value
                if (firstQuals[firstElement.getOffset()] < 0
                        || firstQuals[firstElement.getOffset()] > QualityUtils.MAX_SAM_QUAL_SCORE) {
                    firstQuals[firstElement.getOffset()] = QualityUtils.MAX_SAM_QUAL_SCORE;
                }
                secondQuals[secondElement.getOffset()] = 0;
            } else {
                // if not, we lost confidence in the one with higher quality
                if (firstElement.getQual() >= secondElement.getQual()) {
                    firstQuals[firstElement.getOffset()] =
                            (byte) (SAMTOOLS_OVERLAP_LOW_CONFIDENCE * firstQuals[firstElement.getOffset()]);
                    secondQuals[secondElement.getOffset()] = 0;
                } else {
                    secondQuals[secondElement.getOffset()] =
                            (byte) (SAMTOOLS_OVERLAP_LOW_CONFIDENCE * secondQuals[secondElement.getOffset()]);
                    firstQuals[firstElement.getOffset()] = 0;
                }
            }
            firstElement.getRead().setBaseQualities(firstQuals);
            secondElement.getRead().setBaseQualities(secondQuals);
        }
    }
}
