package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs;

import org.apache.commons.lang3.tuple.MutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.util.*;

/**
 * K-best sub-haplotype finder that selects the best solutions out of a collection of sub-haplotype finders.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public class AggregatedSubHaplotypeFinder<F extends KBestSubHaplotypeFinder> implements KBestSubHaplotypeFinder {

    /**
     * Collection of subFinders that provided the actual solutions.
     */
    final Collection<F> subFinders;

    /**
     *  Flag indicating whether the sub-finders have been processed or not.
     */
    private boolean processedSubFinders = false;

    /**
     * Holds the number of k-best solution that this finder would ever return.
     */
    private int count = 0;

    /**
     * Holds the best {@code i} paths to the sink so far calculated where {@code i+1} is the length of this list.
     *
     * <p>As more results are requested the array will grow. All positions and solutions are
     * calculated up to {@code i}</p>.
     */
    private List<KBestHaplotype> rankedSubHaplotype;

    /**
     * Priority queue with next best haplotype solution from each sub-finder; previous ones are
     * already part {@link #rankedSubHaplotype}.
     */
    private PriorityQueue<MyKBestHaplotypeResult> nextBestSubHaplotypes;

    /**
     * Creates a new aggregated sub-haplotype finder given its sub-finders.
     * @param finders set of sub-finders.
     */
    public AggregatedSubHaplotypeFinder(final Collection<F> finders) {
        Utils.nonNull(finders, "finder collection cannot be null");
        this.subFinders = finders;
    }

    @Override
    public String id() {
        final StringBuilder resultBuilder = new StringBuilder();
        for (final KBestSubHaplotypeFinder subFinder : subFinders) {
            resultBuilder.append(subFinder.id());
        }
        return resultBuilder.toString();
    }

    @Override
    public String label() {
        return "&lt;OR&gt;";
    }

    @Override
    public Set<Pair<? extends KBestSubHaplotypeFinder, String>> subFinderLabels() {
        final int subFinderCount = subFinders.size();
        final String edgeCost = String.format("%.2f", -Math.log10(subFinderCount));
        final Set<Pair<? extends  KBestSubHaplotypeFinder, String>> result = new LinkedHashSet<>(subFinderCount);
        for (final KBestSubHaplotypeFinder subFinder : subFinders) {
            result.add(new MutablePair<>(subFinder, edgeCost));
        }
        return result;
    }

    @Override
    public int getCount() {
        processSubFindersIfNeeded();
        return count;
    }

    @Override
    public double score(final byte[] bases, final int offset, final int length) {
        Utils.nonNull(bases, "bases cannot be null");
        ParamUtils.isPositiveOrZero(offset, "the offset cannot be negative");
        ParamUtils.isPositiveOrZero(length, "the length cannot be negative");
        Utils.validateArg(offset + length <= bases.length, "the offset and length go beyond the array size");
        for (final KBestSubHaplotypeFinder subFinder : subFinders) {
            final double score = subFinder.score(bases,offset,length);
            if (!Double.isNaN(score)) {
                return score;
            }
        }
        return Double.NaN;
    }

    private void processSubFindersIfNeeded() {
        if (processedSubFinders) {
            return;
        }

        long count = 0;
        nextBestSubHaplotypes = new PriorityQueue<>(subFinders.size());
        for (final KBestSubHaplotypeFinder finder : subFinders) {
            final int finderCount = finder.getCount();
            if (finderCount != 0) {
                count += finderCount;
                nextBestSubHaplotypes.add(new MyKBestHaplotypeResult(finder,0));
            }
        }

        this.count = (int) Math.min(Integer.MAX_VALUE, count);

        rankedSubHaplotype = new ArrayList<>(10);
        processedSubFinders = true;
    }

    @Override
    public KBestHaplotype getKBest(final int k) {
        ParamUtils.isPositiveOrZero(k, "k cannot be negative");
        processSubFindersIfNeeded();
        Utils.validateArg(k < count, "k must be less than the count");
        if (k < rankedSubHaplotype.size()) {
            return rankedSubHaplotype.get(k);
        }

        for (int i = rankedSubHaplotype.size(); i <= k; i++) {
            // since k < possibleHaplotypeCount is guarantee no to be empty.
            if (nextBestSubHaplotypes.isEmpty()) {
                throw new IllegalStateException("what the heck " + k + ' ' + count);
            }
            final MyKBestHaplotypeResult nextResult = nextBestSubHaplotypes.remove();
            nextResult.rank = i;
            rankedSubHaplotype.add(nextResult);
            final int subRank = nextResult.result.rank();

            // if there is no further solution from the same child we cannot add another solution from that child.
            if (subRank + 1 >= nextResult.subFinder.getCount()) {
                continue;
            }
            nextBestSubHaplotypes.add(new MyKBestHaplotypeResult(nextResult.subFinder, subRank + 1));
        }
        return rankedSubHaplotype.get(k);
    }

    @Override
    public boolean isReference() {
        return false;
    }

    /**
     * Custom implementation of {@link KBestHaplotype} to encapsulate sub-finder results.
     */
    private static final class MyKBestHaplotypeResult extends KBestHaplotype {

        private final KBestSubHaplotypeFinder subFinder;

        private final KBestHaplotype result;

        private int rank;

        private MyKBestHaplotypeResult(final KBestSubHaplotypeFinder finder, final int rank) {
            this.subFinder = finder;
            this.result = finder.getKBest(rank);
            this.rank = -1;
        }

        @Override
        public SeqGraph graph() {
            return result.graph();
        }

        @Override
        public double score() {
            return result.score();
        }

        @Override
        public boolean isReference() {
            return result.isReference();
        }

        @Override
        public int rank() {
            return rank;
        }

        @Override
        protected SeqVertex head() {
            return result.head();
        }

        @Override
        protected KBestHaplotype tail() {
            return result.tail();
        }
    }
}
