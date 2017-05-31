package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import com.google.common.annotations.VisibleForTesting;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;

/**
 * generic utility class that counts kmers
 *
 * Basically you add kmers to the counter, and it tells you how many occurrences of each kmer it's seen.
 */
public final class KMerCounter {

    /**
     * A map of for each kmer to its num occurrences in addKmers
     */
    private final Map<Kmer, CountedKmer> countsByKMer = new HashMap<>();
    private final int kmerLength;

    /**
     * Create a new kmer counter
     *
     * @param kmerLength the length of kmers we'll be counting to error correct, must be >= 1
     */
    public KMerCounter(final int kmerLength) {
        Utils.validateArg( kmerLength > 0, () -> "kmerLength must be > 0 but got " + kmerLength);
        this.kmerLength = kmerLength;
    }

    /**
     * Get the count of kmer in this kmer counter
     * @param kmer a non-null counter to get
     * @return a positive integer
     */
    public int getKmerCount(final Kmer kmer) {
        Utils.nonNull(kmer, "kmer cannot be null");
        final CountedKmer counted = countsByKMer.get(kmer);
        return counted == null ? 0 : counted.count;
    }

    /**
     * Get an unordered collection of the counted kmers in this counter
     * @return a non-null collection
     */
    public Collection<CountedKmer> getCountedKmers() {
        return countsByKMer.values();
    }

    /**
     * Get kmers that have minCount or greater in this counter
     * @param minCount only return kmers with count >= this value
     * @return a non-null collection of kmers
     */
    public Collection<Kmer> getKmersWithCountsAtLeast(final int minCount) {
        final List<Kmer> result = new LinkedList<>();
        for ( final CountedKmer countedKmer : getCountedKmers() ) {
            if ( countedKmer.count >= minCount ) {
                result.add(countedKmer.kmer);
            }
        }
        return result;
    }

    /**
     * Remove all current counts, resetting the counter to an empty state
     */
    public void clear() {
        countsByKMer.clear();
    }

    /**
     * Add a kmer that occurred kmerCount times
     *
     * @param kmer a kmer
     * @param kmerCount the number of occurrences
     */
    public void addKmer(final Kmer kmer, final int kmerCount) {
        Utils.validateArg(kmer.length() == kmerLength, () -> "bad kmer length " + kmer + " expected size " + kmerLength);
        Utils.validateArg( kmerCount >= 0, () -> "bad kmerCount " + kmerCount);

        CountedKmer countFromMap = countsByKMer.get(kmer);
        if ( countFromMap == null ) {
            countFromMap = new CountedKmer(kmer);
            countsByKMer.put(kmer, countFromMap);
        }
        countFromMap.count += kmerCount;
    }

    @Override
    public String toString() {
        final StringBuilder b = new StringBuilder("KMerCounter{");
        b.append("counting ").append(countsByKMer.size()).append(" distinct kmers");
        b.append("\n}");
        return b.toString();
    }

    static final class CountedKmer implements Comparable<CountedKmer> {
        final Kmer kmer;
        int count = 0;

        private CountedKmer(final Kmer kmer) {
            this.kmer = kmer;
        }

        public Kmer getKmer() {
            return kmer;
        }

        public int getCount() {
            return count;
        }

        @Override
        public String toString() {
            return "CountedKmer{" +
                    "kmer='" + kmer + '\'' +
                    ", count=" + count +
                    '}';
        }

        @Override
        public int compareTo(final CountedKmer o) {
            return Integer.compare(o.count, count);
        }
    }

    @VisibleForTesting
    void addKmer(final String rawKmer, final int kmerCount) {
        addKmer(new Kmer(rawKmer), kmerCount);
    }

    @VisibleForTesting
    void addKmers(final String... kmers) {
        for ( final String kmer : kmers ) {
            addKmer(kmer, 1);
        }
    }
}
