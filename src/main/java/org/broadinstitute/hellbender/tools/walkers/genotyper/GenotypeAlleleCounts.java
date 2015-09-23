package org.broadinstitute.hellbender.tools.walkers.genotyper;

import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.utils.MathUtils;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * Collection of allele counts for a genotype. It encompasses what alleles are present in the genotype and in what number.</p>
 *
 * <p>Alleles are represented herein by their indices running from <b>0</b> to <b>N-1</b> where <i>N</i> is the number of alleles.</p>
 *
 * <p>Each allele present in a genotype (count != 0) has a <i>rank</i>, that is the 0-based ordinal of
 * that allele amongst the ones present in the genotype as sorted by their index.</p>
 *
 * <p>For example:</p>
 *
 * <p><b>0/0/2/2</b> has two alleles with indices <b>0</b> and <b>2</b>, both with count 2.
 * The rank of <b>0</b> is <i>0</i> whereas the rank of <b>2</b> is <i>1</i>.</p>
 *
 * <p><b>2/4/4/7</b> has three alleles with indices <b>2</b>, <b>4</b> and <b>7</b>. <b>2</b> and <b>7</b> have count 1 whereas <b>4</b> has count 2.
 * The rank of <b>2</b> is <i>0</i>, the rank of <b>4</b> is <i>1</i>. and the rank of <b>7</b> is <i>2</i>.</p>
 *
 * <p>In contrast, in both examples above both <b>3</b> and <b>10</b> (and many others) are absent thus they have no rank (represented by <i>-1</i> whenever applies).</p>
 *
 * <p>{@link GenotypeAlleleCounts} instances have themselves their own index (returned by {@link #index() index()}, that indicate their 0-based ordinal within the possible genotype combinations with the same ploidy.</p>
 *
 * <p>For example, for ploidy 3:</p>
 *
 * <table>
 *     <th>Index</th><th>Genotype</th>
 *     <tr><td>0</td><td><b>0/0/0</b></td></tr>
 *     <tr><td>1</td><td><b>0/0/1</b></td></tr>
 *     <tr><td>2</td><td><b>0/1/1</b></td></tr>
 *     <tr><td>3</td><td><b>1/1/1</b></td></tr>
 *     <tr><td>4</td><td><b>0/0/2</b></td></tr>
 *     <tr><td>6</td><td><b>0/1/2</b></td></tr>
 *     <tr><td>7</td><td><b>1/1/2</b></td></tr>
 *     <tr><td>8</td><td><b>0/2/2</b></td></tr>
 *     <tr><td>9</td><td><b>1/2/2</b></td></tr>
 *     <tr><td>10</td><td><b>2/2/2</b></td></tr>
 *     <tr><td>11</td><td><b>0/0/3</b></td></tr>
 *     <tr><td>12</td><td><b>0/1/3</b></td></tr>
 *     <tr><td>13</td><td><b>1/1/3</b></td></tr>
 *     <tr><td>14</td><td><b>0/2/3</b></td></tr>
 *     <tr><td>15</td><td><b>1/2/3</b></td></tr>
 *     <tr><td>16</td><td><b>2/2/3</b></td></tr>
 *     <tr><td>17</td><td><b>0/3/3</b></td></tr>
 *     <tr><td>...</td><td>...</td></tr>
 * </table>
 *
 * The total number of possible genotypes is only bounded by the maximum allele index.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class GenotypeAlleleCounts implements Comparable<GenotypeAlleleCounts> {

    private double log10CombinationCount;

    /**
     * The ploidy of the genotype.
     */
    private final int ploidy;

    /**
     * Sorted array of integer pairs as described in {@link #GenotypeAlleleCounts(int, int, int...)}.
     */
    private int[] sortedAlleleCounts;

    /**
     * Number of different alleles in the genotype.
     */
    private int distinctAlleleCount;

    /**
     * Index of this genotype within genotypes of the same ploidy.
     */
    private int index;

    /**
     * Creates a new unphased genotype.
     *
     * <p>This method assumes that the invoker is passing a well formatted and sorted allele frequency array.
     * Not checks are done for the sake of performance.</p>
     *
     * <p>
     *     The input argument {@code sortedAlleleCounts} list the index of alleles included in the unphased genotype
     *     and their frequency in the genotype in a single array using consecutive pairs:<br/>
     *
     *     <pre> [allele_1,freq_1,allele_2,freq_2, ... , allele_i, freq_i, ... , allele_n, freq_n]</pre>
     *
     *     <br/>
     *     No entry can have frequency == 0 (these must be omitted) and entries are sorted by allele index without
     *     any repetitions so that if <i>i < j</i> then <i>allele_i < allele_j</i>.
     *
     * </p>
     *
     * <p>
     *     The {@code ploidy} provided must be equal to the sum of all frequencies in {@code sortedAlleleCounts}
     * </p>
     * @param ploidy the genotype ploidy.
     * @param sortedAlleleCounts sorted allele counts following the restrictions above.
     * @param index the genotype index.
     */
    private GenotypeAlleleCounts(final int ploidy, final int index, final int... sortedAlleleCounts) {
        this(ploidy, index, sortedAlleleCounts, sortedAlleleCounts.length >> 1, -1);
    }

    private GenotypeAlleleCounts(final int ploidy, final int index, final int[] sortedAlleleCounts, final int distinctAlleleCount, final double log10CombinationCount){
        this.ploidy = ploidy;
        this.index = index;
        this.sortedAlleleCounts = sortedAlleleCounts;
        this.distinctAlleleCount = distinctAlleleCount;
        this.log10CombinationCount = log10CombinationCount;
    }

    /**
     * Returns the log10 of the number of possible allele combinations that would give raise to this allele count.
     * @return 0 or less.
     */
    public double log10CombinationCount() {
        if (log10CombinationCount == -1) {
            return log10CombinationCount = calculateLog10CombinationCount();
        } else {
            return log10CombinationCount;
        }
    }

    /**
     * Calculates log10 combination count.
     *
     * @return 0 or less.
     */
    private double calculateLog10CombinationCount() {
        if (ploidy <= 1) {
            return 0;
        } else {
            final int[] counts = new int[distinctAlleleCount];
            for (int i = 0, j = 1; i < distinctAlleleCount; i++, j+=2) {
                counts[i] = sortedAlleleCounts[j];
            }
            return MathUtils.log10MultinomialCoefficient(ploidy, counts);
        }
    }

    /**
     * Returns the genotype's ploidy.
     * @return 0 or greater.
     */
    public int ploidy() {
        return ploidy;
    }

    /**
     * Increases the allele counts a number of times.
     *
     * <p>
     *     This method must not be invoked on cached genotype-allele-counts that are meant to remain constant,
     *     such as the ones contained in {@link GenotypeLikelihoodCalculators#genotypeTableByPloidy}.
     * </p>
     *
     * @param times the number of times to increase.
     *
     * @throws IllegalArgumentException if {@code times} is negative.
     */
    protected void increase(final int times) {
        for (int i = 0; i < times; i++) {
            increase();
        }
    }

    /**
     * Updates the genotype counts to match the next genotype.
     *
     * <p>
     *     This method must not be invoked on cached genotype-allele-counts that are meant to remain constant,
     *     such as the ones contained in {@link GenotypeLikelihoodCalculators#genotypeTableByPloidy}
     * </p>
     */
    protected void increase() {
        // if the ploidy is zero there is only one possible genotype.
        if (distinctAlleleCount == 0) {
            return;
        }

        // Worth make this case faster.
        if (distinctAlleleCount == 1) {
            if (ploidy == 1) {
                sortedAlleleCounts[0]++;
            } else {
                if (sortedAlleleCounts.length < 4) {
                    sortedAlleleCounts = Arrays.copyOf(sortedAlleleCounts, 4);
                }
                sortedAlleleCounts[2] = sortedAlleleCounts[0] + 1;
                sortedAlleleCounts[3] = 1;
                sortedAlleleCounts[0] = 0;
                sortedAlleleCounts[1] = ploidy - 1;
                distinctAlleleCount = 2;
            }
        } else {
            // Now, all the following ifs are just the way to avoid working with dynamically sizing List<int[]>
            // as the final size of the resulting new sorted-allele-counts array varies depending on the situation.
            // this is considerably faster and the logic complexity would not be that different actually so it is worth
            // the if indentations.
            //
            // Notice that at this point distinctAlleleCount >= 2 thus sortedAlleleCounts.length >= 4.
            //
            // We only need to look at the two lowest allele indices to decide what to do.

            final int allele0 = sortedAlleleCounts[0];
            final int freq0 = sortedAlleleCounts[1];
            final int allele1 = sortedAlleleCounts[2];
            final int allele0Plus1 = allele0 + 1;
            final boolean allele0And1AreConsecutive = allele0Plus1 == allele1;
            final int[] newSortedAlleleCounts;
            // The rest of the sorted allele counts array contains junk
            final int sortedAlleleCountsLength = distinctAlleleCount << 1;


            if (freq0 == 1) {   // in this case allele0 wont be present in the result and all is frequency should go to allele0 + 1.
                if (allele0And1AreConsecutive) {  // need just to remove the first allele and add 1 to the frequency of the second (freq1 += 1).
                    System.arraycopy(sortedAlleleCounts, 2, sortedAlleleCounts, 0, sortedAlleleCountsLength - 2); // shift left the first component away.
                    sortedAlleleCounts[1]++; // freq1 has become freq0.
                    distinctAlleleCount--;
                } else  // just need to mutate allele0 to allele0 + 1.
                {
                    sortedAlleleCounts[0] = allele0Plus1;
                }
            } else { // && freq0 > 1 as per sortedAlleleCounts format restrictions. In this case allele0 will mutated to '0' with frequency decreased by 1.
                if (allele0And1AreConsecutive) { // we don't need to add a component for allele0 + 1 since it already exists.
                    sortedAlleleCounts[0] = 0;
                    sortedAlleleCounts[1] = freq0 - 1;
                    sortedAlleleCounts[3]++;
                } else { // we need to insert allele0 + 1 in the sorted-allele-counts array and give it frequency 1.
                    if (sortedAlleleCounts.length < sortedAlleleCountsLength + 2) // make room for the new component.
                    {
                        sortedAlleleCounts = Arrays.copyOf(sortedAlleleCounts, sortedAlleleCountsLength + 2);
                    }
                    System.arraycopy(sortedAlleleCounts, 2, sortedAlleleCounts, 4, sortedAlleleCountsLength - 2);
                    sortedAlleleCounts[0] = 0;
                    sortedAlleleCounts[1] = freq0 - 1;
                    sortedAlleleCounts[2] = allele0Plus1;
                    sortedAlleleCounts[3] = 1;
                    distinctAlleleCount++;
                }
            }
        }
        index++;
        log10CombinationCount = -1;
    }

    /**
     * Calculates the next genotype in likelihood indexing order.
     * @return never null.
     */
    protected GenotypeAlleleCounts next() {
        // if the ploidy is zero there is only one possible genotype.
        if (distinctAlleleCount == 0) {
            return this;
        }

        // Worth make this case faster.
        if (distinctAlleleCount == 1) {
            if (ploidy == 1)  // A -> B , D -> E etc...
            {
                return new GenotypeAlleleCounts(1, index + 1, sortedAlleleCounts[0] + 1, 1);
            } else  // AAAAA -> AAAAB, DDD -> AAE etc...
            {
                return new GenotypeAlleleCounts(ploidy, index + 1, 0, ploidy - 1, sortedAlleleCounts[0] + 1, 1);
            }
        }
        // Now, all the following ifs are just the way to avoid working with dynamically sizing List<int[]>
        // as the final size of the resulting new sorted-allele-counts array varies depending on the situation.
        // this is considerably faster and the logic complexity would not be that different actually so it is worth
        // the if indentations.
        //
        // Notice that at this point distinctAlleleCount >= 2 thus sortedAlleleCounts.length >= 4.
        //
        // We only need to look at the two lowest allele indices to decide what to do.

        final int allele0 = sortedAlleleCounts[0];
        final int freq0 = sortedAlleleCounts[1];
        final int allele1 = sortedAlleleCounts[2];
        final int allele0Plus1 = allele0 + 1;
        final boolean allele0And1AreConsecutive = allele0Plus1 == allele1;
        final int[] newSortedAlleleCounts;
        // The rest of the sorted allele counts array contains junk
        final int sortedAlleleCountsLength = distinctAlleleCount << 1;

        if (freq0 == 1) {   // in this case allele0 wont be present in the result and all is frequency should go to allele0 + 1.
            if (allele0And1AreConsecutive) {  // need just to remove the first allele and 1 to the frequency of the second (freq1 += 1).
                newSortedAlleleCounts = Arrays.copyOfRange(sortedAlleleCounts, 2, sortedAlleleCountsLength);
                newSortedAlleleCounts[1]++;
            } else {  // just need to mutate allele0 to allele0 + 1.
                newSortedAlleleCounts = Arrays.copyOf(sortedAlleleCounts, sortedAlleleCountsLength);
                newSortedAlleleCounts[0] = allele0Plus1;
                // newSortedAlleleCounts[1] = 1; // :) no need to do it because it is already the case (freq0 == 1).
            }
        } else { // && freq0 > 1 as per sortedAlleleCounts format restrictions. In this case allele0 will muttated to '0' with frequency decreased by 1.
            if (allele0And1AreConsecutive) { // we don't need to add a component for allele0 + 1 since it already exists.
                newSortedAlleleCounts = sortedAlleleCounts.clone();
                newSortedAlleleCounts[0] = 0;
                newSortedAlleleCounts[1] = freq0 - 1;
                newSortedAlleleCounts[3]++;
            } else { // we need to insert allele0 + 1 in the sorted-allele-counts array.
                newSortedAlleleCounts = new int[sortedAlleleCountsLength + 2];
                newSortedAlleleCounts[0] = 0;
                newSortedAlleleCounts[1] = freq0 - 1;
                newSortedAlleleCounts[2] = allele0Plus1;
                newSortedAlleleCounts[3]++; // = 1 as the array was freshly created with 0s.
                System.arraycopy(sortedAlleleCounts, 2, newSortedAlleleCounts, 4, sortedAlleleCountsLength - 2);
            }
        }
        return new GenotypeAlleleCounts(ploidy, index + 1, newSortedAlleleCounts);
    }

    /**
     * Returns the number of different alleles that participate in the genotype.
     *
     * @return 0 or greater.
     */
    public int distinctAlleleCount() {
        return distinctAlleleCount;
    }

    /**
     * Returns the index of the allele from its rank in the genotype.
     *
     * @param rank the query rank.
     *
     * @throws IllegalArgumentException if the {@code rank} provided is outside the valid range [0,{@link #distinctAlleleCount()}).
     *
     * @return 0 or greater.
     */
    public int alleleIndexAt(final int rank) {
        if (rank < 0 || rank >= distinctAlleleCount) {
            throw new IllegalArgumentException("the requested rank " + rank + " is out of range [0," + distinctAlleleCount + ")");
        }
        return sortedAlleleCounts[rank << 1];
    }

    /**
     * Returns the rank of an allele in the genotype by its index.
     *
     * @param index the target index.
     *
     * @throws IllegalArgumentException if {@code index} is less that 0. Indices can be arbitrarily large.
     *
     * @return -1 or less if the allele index is not present in the genotype, 0 to {@link #distinctAlleleCount()} - 1 otherwise.
     *   If negative, the absolute value can be used to determine where would be that index inserted within {@code [0,{@link #distinctAlleleCount()}]} as
     *   {@code - result - 1}.
     *
     */
    public int alleleRankFor(final int index) {
        if (index < 0) {
            throw new IllegalArgumentException("the index must be 0 or greater");
        }
        return alleleIndexToRank(index, 0, distinctAlleleCount);
    }

    /**
     * Generates a string that would represent the unphased genotype with this allele counts.
     *
     * <p>
     *     In this string allele calls appear in alleleIndex order with as many repeats as copies of each allele. So
     *     for example:<br/>
     *     <pre>
     *         0         # haploid reference.
     *         0/0       # typical diploid calls
     *         0/1
     *         1/1
     *         0/0/1/3/3 # pentaploid with to ref, one first alt. and 2 third alt. allele
     *     </pre>
     *
     * </p>
     *
     * @return never {@code null}.
     */
    public String toUnphasedGenotypeString() {
        if (ploidy == 0) {
            return "";
        }
        final StringBuilder sb = new StringBuilder(distinctAlleleCount * 3);
        for (int i = 0; i < distinctAlleleCount; i += 2) {
            final int alleleIndex = sortedAlleleCounts[i];
            final int alleleCount = sortedAlleleCounts[i + 1];
            for (int j = 0; j < alleleCount; j++) {
                sb.append(alleleIndex).append('/');
            }

        }
        sb.setLength(sb.length() - 1);
        return sb.toString();
    }

    @Override
    public String toString() {
        // Perhaps we should change in the future, but the unphased genotype representation seems to be
        // a good one.
        return toUnphasedGenotypeString();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public boolean equals(final Object o) {
        if (o instanceof GenotypeAlleleCounts) {
            return equals((GenotypeAlleleCounts) o);
        } else {
            return false;
        }
    }

    /**
     * Compares with another genotype.
     * @param o the other genotype.
     * @return never {@code null}.
     */
    public boolean equals(final GenotypeAlleleCounts o) {
        if (o == this) {
            return true;
        }
        if (o == null) {
            return false;
        }
        if (ploidy != o.ploidy) {
            return false;
        }
        return Arrays.equals(sortedAlleleCounts, o.sortedAlleleCounts);
    }

    /**
     * Returns the index of this genotype allele count within all possible genotypes with the same ploidy.
     *
     * @return 0 or greater.
     */
    public int index() {
        return index;
    }

    /**
     * Compares to genotypes.
     *
     * <p>A genotype with larger ploidy is considered greater than one with a lower ploidy. If both genotypes have
     * the same ploidy, then the genotype with the largest allele index or largest count if these are the same</p>.
     *
     * @param other genotype to compare to.
     *
     * @throws IllegalArgumentException if {@code other} is {@code null}.
     *
     * @return 0 if both genotypes are equivalent, < 0 if this genotype is less than {@code other} and > 0
     * if this genotype is greater than {@code other}.
     */
    @Override
    public int compareTo(final GenotypeAlleleCounts other) {
        if (other == this) {
            return 0;
        }
        if (other == null) {
            throw new IllegalArgumentException("input genotype cannot be null");
        }
        if (other.ploidy == ploidy) {
            return index - other.index;
        } else {
            return ploidy - other.ploidy;
        }
    }

    @Override
    public int hashCode() {
        return ((31 + ploidy) * 31 ) + index;
    }

    /**
     * Implements binary search across allele indexes.
     * @param index the target index.
     * @param from first inclusive possible rank.
     * @param to last exclusive possible rank.
     * @return -1 or less if the allele index is not in the genotype false otherwise. You can obtain
     *  the potential insertion point (within the interval [from,to]) as {@code -result - 1}
     */
    private int alleleIndexToRank(final int index,final int from, final int to) {
        if (to <= from) {
            return -from - 1;
        }
        if (from == to - 1) {
            final int onlyIndex = sortedAlleleCounts[from << 1];
            return onlyIndex == index ? from : (onlyIndex > index) ? -from - 1 : -to - 1;
        }

        final int mid = (to + from) >> 1;
        final int midIndex = sortedAlleleCounts[mid << 1];
        if (midIndex == index) {
            return mid;
        } else if (midIndex < index) {
            return alleleIndexToRank(index, mid + 1, to);
        } else {
            return alleleIndexToRank(index, 0, mid);
        }
    }

    /**
     * Returns the count of an allele in the genotype given is rank in the genotype (not the allele index itself).
     *
     * @param rank of the requested allele within the genotype.
     *
     * @throws IllegalArgumentException if {@code rank} is out the the valid range [0,{@link #distinctAlleleCount})
     *
     * @return 1 or greater.
     */
    public int alleleCountAt(final int rank) {
        if (rank < 0 || rank >= distinctAlleleCount) {
            throw new IllegalArgumentException("the rank is out of range");
        }
        return sortedAlleleCounts[(rank << 1) + 1];
    }

    /**
     * Checks whether this genotype contain at least one call on a particular allele index.
     *
     * @param index the target allele.
     *
     * @throws IllegalArgumentException if {@code index} is negative.
     *
     * @return {@code true} iff the genotype contains that allele index.
     */
    public boolean containsAllele(final int index) {
        return alleleRankFor(index) >= 0;
    }

    /**
     * Returns the count of an allele in the genotype given it index.
     *
     * @return 0 if the allele is not present in the genotype, 1 or more otherwise.
     */
    public int alleleCountFor(final int index) {
        final int rank = alleleRankFor(index);
        return rank < 0 ? 0 : alleleCountAt(rank);
    }

    /**
     * Returns the allele counts for each allele index to maximum.
     * @param maximumAlleleIndex the maximum allele index required.
     * @throws IllegalArgumentException if {@code maximumAlleleIndex} is less than 0.
     * @return never {@code null}, an array of exactly {@code maximumAlleleIndex + 1} positions with the counts
     * of each allele where the position in the array is equal to its index.
     */
    public int[] alleleCountsByIndex(final int maximumAlleleIndex) {
        if (maximumAlleleIndex < 0) {
            throw new IllegalArgumentException("the requested allele count cannot be less than 0");
        }
        final int[] result = new int[maximumAlleleIndex + 1];
        copyAlleleCountsByIndex(result, 0, 0, maximumAlleleIndex);
        return result;
    }


    private void copyAlleleCountsByIndex(final int[] dest, final int offset, final int minimumAlleleIndex, final int maximumAlleleIndex) {

        // First we determine what section of the sortedAlleleCounts array contains the counts of interest,
        // By the present allele rank range of interest.
        final int minimumAlleleRank = alleleRankFor(minimumAlleleIndex);
        final int maximumAlleleRank = alleleRankFor(maximumAlleleIndex);

        // If the min or max allele index are absent (returned rank < 0) we note where the would be inserted; that
        // way we avoid going through the rest of positions in the sortedAlleleCounts array.
        // The range of interest is then [startRank,endRank].
        final int startRank = minimumAlleleRank < 0 ? - minimumAlleleRank - 1 : minimumAlleleRank;
        final int endRank = maximumAlleleRank < 0 ? - maximumAlleleRank - 2 : maximumAlleleRank;

        // Iteration variables:
        int nextIndex = minimumAlleleIndex; // next index that we want to output the count for.
        int nextRank = startRank; // next rank to query in sortedAlleleCounts.
        int nextSortedAlleleCountsOffset = nextRank << 1; // offset in sortedAlleleCounts where the info is present for the next rank.
        int nextDestOffset = offset; // next offset in destination array where to set the count for the nextIndex.

        while (nextRank++ <= endRank) {
            final int alleleIndex = sortedAlleleCounts[nextSortedAlleleCountsOffset++];
            // fill non-present allele counts with 0s.
            while (alleleIndex > nextIndex) {
                dest[nextDestOffset++] = 0;
                nextIndex++;
            }
            // It is guaranteed that at this point alleleIndex == nextIndex
            // thanks to the condition of the enclosing while: there must be at least one index of interest that
            // is present in the remaining (nextRank,endRank] interval as otherwise endRank would be less than nextRank.
            dest[nextDestOffset++] = sortedAlleleCounts[nextSortedAlleleCountsOffset++];
            nextIndex++;
        }
        // Finally we take care of trailing requested allele indices.
        while (nextIndex++ <= maximumAlleleIndex) {
            dest[nextDestOffset++] = 0;
        }
    }

    /**
     * Copies the sorted allele counts into an array.
     *
     * <p>
     *     Sorted allele counts are disposed as an even-sized array where even positions indicate the allele index and
     *     the following odd positions the number of copies of that allele in this genotype allele count:
     * </p>
     * <p><pre>
     *     [ allele_0, freq_0, allele_1, freq_1 ... ]
     * </pre></p>
     *
     * <p>
     *     With {@code offset} you can indicate an alternative first position in the destination array.
     * </p>
     *
     * @param dest where to copy the counts.
     * @param offset starting position.
     *
     * @throws IllegalArgumentException if {@code dest} is {@code null}, {@code offset} is less than 0
     *   or {@code dest} is not large enough considering the number of alleles present in this genotype
     *   allele counts and the {@code offset} provided. A total of
     *   <code>{@link #distinctAlleleCount()} * 2 positions</code>
     *   are required for the job.
     */
    public void copyAlleleCounts(final int[] dest, final int offset) {
        if (dest == null) {
            throw new IllegalArgumentException("the destination cannot be null");
        }
        if (offset < 0) {
            throw new IllegalArgumentException("the offset cannot be negative");
        }
        final int sortedAlleleCountsLength = distinctAlleleCount << 1;
        if (offset + sortedAlleleCountsLength > dest.length) {
            throw new IllegalArgumentException("the input array does not have enough capacity");
        }
        System.arraycopy(sortedAlleleCounts, 0, dest, offset, sortedAlleleCountsLength);
    }

    /**
     * Instantiates the first genotype possible provided a total ploidy.
     * @param ploidy the ploidy of the genotype.
     *
     * @throws java.lang.IllegalArgumentException if ploidy is less than 0.
     *
     * @return never {@code null}.
     */
    protected static GenotypeAlleleCounts first(final int ploidy) {
        if (ploidy < 0) {
            throw new IllegalArgumentException("the ploidy must be 0 or greater");
        } else if (ploidy == 0) {
            return new GenotypeAlleleCounts(0, 0);
        } else {
            return new GenotypeAlleleCounts(ploidy, 0, 0, ploidy);
        }
    }

    /**
     * Makes the next genotype in likelihood indexing order.
     *
     * @param g the original genotype.
     *
     * @throws IllegalArgumentException if {@code g} is {@code null}.
     *
     * @return never {@code null}.
     */
    public static GenotypeAlleleCounts makeNextGenotype(final GenotypeAlleleCounts g) {
        if (g == null) {
            throw new IllegalArgumentException("the next genotype");
        }
        return g.next();
    }

    /**
     * Returns the largest allele index present in the genotype.
     *
     * @return -1 if there is no alleles (ploidy == 0), 0 or greater otherwise.
     */
    public int maximumAlleleIndex() {
        if (distinctAlleleCount == 0) {
            return -1;
        } else {
            return sortedAlleleCounts[(distinctAlleleCount - 1) << 1];
        }
    }

    /**
     * Returns the smallest allele index present in the genotype.
     *
     * @return -1 if there is no allele (ploidy == 0), 0 or greater otherwise.
     */
    public int minimumAlleleIndex() {
        if (distinctAlleleCount == 0) {
            return -1;
        } else {
            return sortedAlleleCounts[0];
        }
    }

    /**
     * Creates an independent copy of this GenotypeAlleleCounts.
     * @return never {@code null}.
     */
    GenotypeAlleleCounts copy() {
        return new GenotypeAlleleCounts(this.ploidy, this.index, this.sortedAlleleCounts.clone(), this.distinctAlleleCount, this.log10CombinationCount);
    }

    /**
     * Composes a list with the alleles.
     *
     * @param allelesToUse alleles to use.
     *
     * @throws IllegalArgumentException if {@code allelesToUse} is {@code null},
     *          or does not contain enough elements to accommodate the maximum allele index in this allele-counts.
     *
     * @return never null, but it might be restricted (unmodifiable or non-expandable).
     */
    @SuppressWarnings("unchecked")
    public <T extends Allele> List<T> asAlleleList(final List<T> allelesToUse) {
        if (allelesToUse == null) {
            throw new IllegalArgumentException("the input allele list cannot be null");
        }
        if (allelesToUse.size() < maximumAlleleIndex()) {
            throw new IllegalArgumentException("the provided alleles to use does not contain an element for the maximum allele ");
        }
        if (distinctAlleleCount == 1 ) {
            if (ploidy == 1) {
                return Collections.singletonList(allelesToUse.get(sortedAlleleCounts[0]));
            } else {
                return Collections.nCopies(ploidy, allelesToUse.get(sortedAlleleCounts[0]));
            }
        } else {
            final Allele[] myAlleles = new Allele[ploidy];
            int nextIndex = 0;
            for (int i = 0, ii = 0; i < distinctAlleleCount; i++) {
                final Allele allele = allelesToUse.get(sortedAlleleCounts[ii++]);
                final int repeats = sortedAlleleCounts[ii++];
                for (int j = 0; j < repeats; j++) {
                    myAlleles[nextIndex++] = allele;
                }
            }
            return (List<T>) Arrays.asList(myAlleles);
        }
    }

    /**
     * Returns an array with the allele indices repeated based on the number of occurrences in the genotype.
     *
     * <p>
     *     indices are sorted from the smallest to the greatest.
     * </p>
     *
     * <p>
     *     If a sufficiently large array is provided as {@code dest}, this is used as the destination. Unnecessary
     *     positions at the back of the array are left untouched.
     * </p>
     *
     * <p>
     *     However if {@code dest} is {@code null} or it does not have enough space, a new array of with length equal to
     *     the ploidy will be used and returned instead.
     * </p>
     *
     * @param dest destination array. Can be {@code null} or not have sufficient positions (ploidy); in that case a new
     *             one is created.
     * @return never {@code null}, {@code dest} if sufficiently large otherwise an array of ploidy positions.
     */
    public int[] toAlleleIndicesArray(final int[] dest) {
        final int[] result = dest == null || dest.length < ploidy ? new int[ploidy] : dest;
        int k = 0;
        for (int i = 0,ii = 0; i < distinctAlleleCount; i++) {
            final int index = sortedAlleleCounts[ii++];
            final int repeats = sortedAlleleCounts[ii++];
            for (int j = 0; j < repeats; j++) {
                result[k++] = index;
            }
        }
        return result;
    }
}
