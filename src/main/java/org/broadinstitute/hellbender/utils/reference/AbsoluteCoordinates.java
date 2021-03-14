package org.broadinstitute.hellbender.utils.reference;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.List;
import java.util.function.IntFunction;
import java.util.function.ToIntFunction;

/**
 * Class to translate back and forth from absolute {@code long-typed} base positions to relative ones (the usual contig, position pairs).
 * <p>
 *     If we were to concatenate the sequence of every chromosome in the order these appear in reference dictionary,
 *     each base in the reference with a absolute position or offset from the beginning of this imaginary super-contig.
 * </p>
 * <p>
 *     Some times it might be handy to use this type of address/coordinate rather than the usual contig (name or index) and position
 *     within the contig.
 * </p>
 * <p>
 *     This class implement such a transformation from absolute coordinates efficiently.
 * </p>
 * <p>
 *     Relative to absolute is rather a trivial matter if you keep track of the accumulative number of bases on the reference
 *     right before the start of each contig.
 * </p>
 * <p> Absolute to relative is a bit trickier. A obvious solution would consist in a binary search amogst the contig for the one that
 *     include the absolute base range that encloses the requested position. The cost of this lookup would be O(log C) where C is the number of
 *     contigs in the reference. In some cases like hg38 there are over 3000+ such contigs due to the alt-contigs and decoy, resulting in around 11 iterations
 *     to find the target's location contig</p>
 * <p>
 *     This class has a couple of accelerations:
 *     <p>
 *         First we check whether the target contig is the last returned, which should be quite often the case when accessing the refernce in
 *         sequence.
 *     </p>
 *     <p>
 *         Then we keep a "percentile" table that contains the contig index that would include that percentile absolute position. The granularity
 *         of such table is linked to the number of contigs in the reference with a O(C) additional memory cost.
 *     </p>
 *     <p>
 *         These two acceleration may effectively reduce the look-up cost to O(1) is most scenarios. The draw back is a more complicate code and the
 *         need to deal with float-point arithmetic for the percentile look-up.
 *     </p>
 * </p>
 */
public final class AbsoluteCoordinates {

    /**
     * Length of each contig.
     */
    private final int[] lengths;

    /**
     * Percentile look up table.
     */
    private final int[] percentiles;

    /**
     * Contains the number of bases before the contig with the ith index (0-based).
     */
    private final long[] accumulative;

    /**
     * Total length of the reference.
     */
    private final long total;

    /**
     * Factor to multiply to an absolute position to get its corresponding percentile.
     */
    private final float percentileFactor;

    /**
     * Function that resolves the contig index given its name.
     */
    private final ToIntFunction<String> contigToIndex;

    /**
     * Function that resolves the contig name given its index.
     */
    private final IntFunction<String> indexToContig;
    private int lastCtg;


    private AbsoluteCoordinates(final int[] lengths, final long[] accumulative,
                                final ToIntFunction<String> contigToIndex, final IntFunction<String> indexToContig) {
        this.lengths = lengths;
        this.accumulative = accumulative;
        this.contigToIndex = contigToIndex;
        this.indexToContig = indexToContig;
        this.total = accumulative[accumulative.length - 1];
        this.percentiles = calculatePercentiles(lengths, accumulative, total);
        this.percentileFactor = (percentiles.length - 1) / (float) this.total;
        this.lastCtg = 0;
    }

    private static int[] calculatePercentiles(final int[] lengths, final long[] accumulative, final long total) {
        final int[] result = new int[(accumulative.length << 1) + 1];
        final float fraction = total / (float)(result.length - 1);
        double fractionAccumulator = 0;
        for (int i = 0, j = 0; i < lengths.length; i++) {
            final long accumulativePlusLength = accumulative[i + 1]; // == accumulative[i] + lengths[i];
            while (fractionAccumulator < accumulativePlusLength && j < result.length - 1) {
                result[j++] = i;
                fractionAccumulator += fraction;
            }
        }
        result[result.length - 1] = lengths.length - 1;
        return result;
    }

    public static AbsoluteCoordinates of(final SAMSequenceDictionary dictionary) {
        final ToIntFunction<String> contigToIndex = dictionary::getSequenceIndex;
        final IntFunction<String> indexToContig = i -> dictionary.getSequence(i).getContig();
        final List<SAMSequenceRecord> sequences = dictionary.getSequences();
        final int numberOfSequences = sequences.size();
        final int[] lengths = new int[numberOfSequences];
        final long[] accumulative = new long[numberOfSequences + 1];
        for (int i = 0; i < numberOfSequences; i++) {
            final SAMSequenceRecord sequence = sequences.get(i);
            lengths[i] = sequence.getSequenceLength();
        }
        long leftSum = lengths[0];
        for (int i = 1; i < numberOfSequences; i++) {
            accumulative[i] = leftSum;
            leftSum += lengths[i];
        }
        accumulative[numberOfSequences] = leftSum;
        return new AbsoluteCoordinates(lengths, accumulative, contigToIndex, indexToContig);
    }

    public long toAbsolute(final String ctgName, final int position) {
        return toAbsolute(contigToIndex.applyAsInt(ctgName), position);
    }

    /**
     * Obtains the absolute coordinate for the start position in an interval.
     * @param simpleInterval the target position in relative coordinates
     * @return 1 or greater.
     * @throws IllegalArgumentException no such contig name or tht conting is too small.
     */
    public long toAbsolute(final SimpleInterval simpleInterval) {
        return toAbsolute(simpleInterval.getContig(), simpleInterval.getStart());
    }


    public long toAbsolute(final int ctgIdx, final int position) {
        if (lengths[ctgIdx] < position) {
            throw new IllegalArgumentException("position outside containg contig");
        }
        lastCtg = ctgIdx;
        return accumulative[ctgIdx] + position;
    }

    public static class Relative {
        public final String contig;
        public final int contigIndex;
        public final int position;

        Relative(final String contig, final int contigIndex, final int position) {
            this.contig = contig;
            this.contigIndex = contigIndex;
            this.position = position;
        }
    }

    @FunctionalInterface
    interface RelativeFactory<E> {
        E create(final String ctgName, final int ctgIdx, final int position);
    }

    public SimpleInterval toSimpleInterval(final long absoluteStart, final int length) {
        return toRelative(absoluteStart, (n, i, p) -> new SimpleInterval(n, p, p + length - 1));
    }

    public Relative toRelative(final long absolute) {
        return toRelative(absolute, Relative::new);
    }

    public <E> E toRelative(final long absolute, final RelativeFactory<E> factory) {
        if (absolute < 1) {
            throw new IllegalArgumentException("absolute cannot be less than 1");
        } else if (absolute > total) {
            throw new IllegalArgumentException("absolute is too large");
        } else if (accumulative[lastCtg] < absolute) {
            if (absolute <= accumulative[lastCtg + 1]) {
               return factory.create(indexToContig.apply(lastCtg), lastCtg, (int) (absolute - accumulative[lastCtg]));
            } else {
                return searchRelative(absolute, factory);
            }
        } else {
            return searchRelative(absolute, factory);
        }
    }

    private <E> E searchRelative(final long target, final RelativeFactory<E> factory) {
        final int percentileIndex = (int) (target * percentileFactor);
        if (percentileIndex >= percentiles.length) {
             throw new IllegalArgumentException("xx " + target + " " + total );
        }
        final int percentileCtg = percentiles[percentileIndex];
        if (percentileIndex < percentiles.length - 1) {
            return searchRelative(target, percentileCtg, percentiles[percentileIndex + 1], factory);
        } else {
            return searchRelative(target, percentileCtg, lengths.length - 1, factory);
        }
    }

    private <E> E searchRelative(final long target, final int minCtgIdx, final int maxCtgIdx, final RelativeFactory<E> factory) {
        int i = minCtgIdx, j = maxCtgIdx;
        while (i < j) {
            final int mid = (i + j) >> 1;
            if (accumulative[mid] >= target) {
                j = mid - 1;
            } else if (accumulative[mid + 1] < target) {
                i = mid + 1;
            } else {
                i = mid;
                break;
            }
        }
        lastCtg = i;
        return factory.create(indexToContig.apply(i), i, (int) (target - accumulative[i]));
    }
}
