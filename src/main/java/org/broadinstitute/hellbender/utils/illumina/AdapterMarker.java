package org.broadinstitute.hellbender.utils.illumina;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.CollectionUtil;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.StringUtil;

import java.util.*;
import java.util.concurrent.atomic.AtomicReference;

/**
 * Store one or more AdapterPairs to use to mark adapter sequence of SAMRecords.  This is a very compute-intensive process, so
 * this class implements two heuristics to reduce computation:
 * - Adapter sequences are truncated, and then any adapter pairs that become identical after truncation are collapsed into a single pair.
 * - After a specified number of reads with adapter sequence has been seen, prune the list of adapter pairs to include only the most
 *   frequently seen adapters.  For a flowcell, there should only be a single adapter pair found.
 *
 * Note that the AdapterPair object returned by all the adapterTrim* methods will not be one of the original AdapterPairs
 * passed to the ctor, but rather will be one of the truncated copies.
 */
public class AdapterMarker {
    public static final int DEFAULT_ADAPTER_LENGTH = 30;
    public static final int DEFAULT_PRUNE_ADAPTER_LIST_AFTER_THIS_MANY_ADAPTERS_SEEN = 100;
    public static final int DEFAULT_NUM_ADAPTERS_TO_KEEP = 1;

    // It is assumed that these are set once during execution, before the class is used to mark any adapters, but this is not enforced.
    private int thresholdForSelectingAdaptersToKeep = DEFAULT_PRUNE_ADAPTER_LIST_AFTER_THIS_MANY_ADAPTERS_SEEN;
    private int numAdaptersToKeep = DEFAULT_NUM_ADAPTERS_TO_KEEP;
    private int minSingleEndMatchBases = ClippingUtil.MIN_MATCH_BASES;
    private int minPairMatchBases = ClippingUtil.MIN_MATCH_PE_BASES;
    private double maxSingleEndErrorRate = ClippingUtil.MAX_ERROR_RATE;
    private double maxPairErrorRate = ClippingUtil.MAX_PE_ERROR_RATE;

    // This is AtomicReference because one thread could be matching adapters while the threshold has been crossed in another
    // thread and the array is being replaced.
    private final AtomicReference<AdapterPair[]> adapters = new AtomicReference<>();

    // All the members below are only accessed within a synchronized block.
    private boolean thresholdReached = false;
    private int numAdaptersSeen = 0;
    private final CollectionUtil.DefaultingMap<AdapterPair, Integer> seenCounts = new CollectionUtil.DefaultingMap<>(0);

    /**
     * Truncates adapters to DEFAULT_ADAPTER_LENGTH
     * @param originalAdapters These should be in order from longest & most likely to shortest & least likely.
     */
    public AdapterMarker(final AdapterPair... originalAdapters) {
        this(DEFAULT_ADAPTER_LENGTH, originalAdapters);
    }

    /**
     * @param adapterLength Truncate adapters to this length.
     * @param originalAdapters These should be in order from longest & most likely to shortest & least likely.
     */
    public AdapterMarker(final int adapterLength, final AdapterPair... originalAdapters) {
        // Truncate each AdapterPair to the given length, and then combine any that end up the same after truncation.
        final ArrayList<TruncatedAdapterPair> truncatedAdapters = new ArrayList<>();
        for (final AdapterPair adapter : originalAdapters) {
            final TruncatedAdapterPair truncatedAdapter = makeTruncatedAdapterPair(adapter, adapterLength);
            final int matchingIndex = truncatedAdapters.indexOf(truncatedAdapter);
            if (matchingIndex == -1) {
                truncatedAdapters.add(truncatedAdapter);
            } else {
                final TruncatedAdapterPair matchingAdapter = truncatedAdapters.get(matchingIndex);
                matchingAdapter.setName(matchingAdapter.getName() + "|" + adapter.getName());
            }
        }
        adapters.set(truncatedAdapters.toArray(new AdapterPair[truncatedAdapters.size()]));
    }

    /**
     * After seeing the thresholdForSelectingAdapters number of adapters, keep up to this many of the original adapters.
     */
    public synchronized AdapterMarker setNumAdaptersToKeep(final int numAdaptersToKeep) {
        if (numAdaptersToKeep <= 0) {
            throw new IllegalArgumentException(String.format("numAdaptersToKeep should be positive: %d", numAdaptersToKeep));
        }
        this.numAdaptersToKeep = numAdaptersToKeep;
        return this;
    }

    /**
     * When this number of adapters have been matched, discard the least-frequently matching ones.
     * @param thresholdForSelectingAdaptersToKeep set to -1 to never discard any adapters.
     */
    public synchronized AdapterMarker setThresholdForSelectingAdaptersToKeep(final int thresholdForSelectingAdaptersToKeep) {
        this.thresholdForSelectingAdaptersToKeep = thresholdForSelectingAdaptersToKeep;
        return this;
    }

    /**
     *
     * @param minSingleEndMatchBases When marking a single-end read, adapter must match at least this many bases.
     */
    public synchronized AdapterMarker setMinSingleEndMatchBases(final int minSingleEndMatchBases) {
        this.minSingleEndMatchBases = minSingleEndMatchBases;
        return this;
    }

    /**
     *
     * @param minPairMatchBases When marking a paired-end read, adapter must match at least this many bases.
     */
    public synchronized AdapterMarker setMinPairMatchBases(final int minPairMatchBases) {
        this.minPairMatchBases = minPairMatchBases;
        return this;
    }

    /**
     * @param maxSingleEndErrorRate For single-end read, no more than this fraction of the bases that align with the adapter can
     *                              mismatch the adapter and still be considered an adapter match.
     */
    public synchronized AdapterMarker setMaxSingleEndErrorRate(final double maxSingleEndErrorRate) {
        this.maxSingleEndErrorRate = maxSingleEndErrorRate;
        return this;
    }

    /**
     * @param maxPairErrorRate For paired-end read, no more than this fraction of the bases that align with the adapter can
     *                         mismatch the adapter and still be considered an adapter match.
     */
    public synchronized AdapterMarker setMaxPairErrorRate(final double maxPairErrorRate) {
        this.maxPairErrorRate = maxPairErrorRate;
        return this;
    }

    public AdapterPair adapterTrimIlluminaSingleRead(final SAMRecord read) {
        return adapterTrimIlluminaSingleRead(read, minSingleEndMatchBases, maxSingleEndErrorRate);
    }

    public AdapterPair adapterTrimIlluminaPairedReads(final SAMRecord read1, final SAMRecord read2) {
        return adapterTrimIlluminaPairedReads(read1, read2, minPairMatchBases, maxPairErrorRate);
    }

    /**
     * Overrides defaults for minMatchBases and maxErrorRate
     */
    public AdapterPair adapterTrimIlluminaSingleRead(final SAMRecord read, final int minMatchBases, final double maxErrorRate) {
        final AdapterPair ret = ClippingUtil.adapterTrimIlluminaSingleRead(read, minMatchBases, maxErrorRate, adapters.get());
        if (ret != null) tallyFoundAdapter(ret);
        return ret;
    }

    /**
     * Overrides defaults for minMatchBases and maxErrorRate
     */
    public AdapterPair adapterTrimIlluminaPairedReads(final SAMRecord read1, final SAMRecord read2,
                                                             final int minMatchBases, final double maxErrorRate) {
        final AdapterPair ret = ClippingUtil.adapterTrimIlluminaPairedReads(read1, read2, minMatchBases, maxErrorRate, adapters.get());
        if (ret != null) tallyFoundAdapter(ret);
        return ret;
    }

    /** For unit testing only */
    AdapterPair[] getAdapters() {
        return adapters.get();
    }

    private TruncatedAdapterPair makeTruncatedAdapterPair(final AdapterPair adapterPair, final int adapterLength) {
        return new TruncatedAdapterPair("truncated " + adapterPair.getName(),
                substringAndRemoveTrailingNs(adapterPair.get3PrimeAdapterInReadOrder(), adapterLength),
                substringAndRemoveTrailingNs(adapterPair.get5PrimeAdapterInReadOrder(), adapterLength));
    }

    /**
     * Truncate to the given length, and in addition truncate any trailing Ns.
     */
    private String substringAndRemoveTrailingNs(final String s, int length) {
        length = Math.min(length, s.length());
        final byte[] bytes = StringUtil.stringToBytes(s);
        while (length > 0 && SequenceUtil.isNoCall(bytes[length - 1])) {
            length--;
        }
        return s.substring(0, length);
    }

    /**
     * Keep track of every time an adapter is found, until it is time to prune the list of adapters.
     */
    private void tallyFoundAdapter(final AdapterPair foundAdapter) {
        // If caller does not want adapter pruning, do nothing.
        if (thresholdForSelectingAdaptersToKeep < 1) return;
        synchronized (this) {
            // Already pruned adapter list, so nothing more to do.
            if (thresholdReached) return;

            // Tally this adapter
            seenCounts.put(foundAdapter, seenCounts.get(foundAdapter) + 1);

            // Keep track of the number of times an adapter has been seen.
            numAdaptersSeen += 1;

            // Reached the threshold for pruning the list.
            if (numAdaptersSeen >= thresholdForSelectingAdaptersToKeep) {

                // Sort adapters by number of times each has been seen.
                final TreeMap<Integer, AdapterPair> sortedAdapters = new TreeMap<>(new Comparator<Integer>() {
                    @Override
                    public int compare(final Integer integer, final Integer integer2) {
                        // Reverse of natural ordering
                        return integer2.compareTo(integer);
                    }
                });
                for (final Map.Entry<AdapterPair, Integer> entry : seenCounts.entrySet()) {
                    sortedAdapters.put(entry.getValue(), entry.getKey());
                }

                // Keep the #numAdaptersToKeep adapters that have been seen the most, plus any ties.
                final ArrayList<AdapterPair> bestAdapters = new ArrayList<>(numAdaptersToKeep);
                int countOfLastAdapter = Integer.MAX_VALUE;
                for (final Map.Entry<Integer, AdapterPair> entry : sortedAdapters.entrySet()) {
                    if (bestAdapters.size() >= numAdaptersToKeep) {
                        if (entry.getKey() == countOfLastAdapter) {
                            bestAdapters.add(entry.getValue());
                        } else {
                            break;
                        }
                    } else {
                        countOfLastAdapter = entry.getKey();
                        bestAdapters.add(entry.getValue());
                    }
                }
                // Replace the existing list with the pruned list.
                thresholdReached = true;
                adapters.set(bestAdapters.toArray(new AdapterPair[bestAdapters.size()]));
            }
        }
    }

    private static class TruncatedAdapterPair implements AdapterPair {
        String name;
        final String fivePrime, threePrime, fivePrimeReadOrder;
        final byte[]  fivePrimeBytes, threePrimeBytes, fivePrimeReadOrderBytes;

        private TruncatedAdapterPair(final String name, final String threePrimeReadOrder, final String fivePrimeReadOrder) {
            this.name = name;
            this.threePrime = threePrimeReadOrder;
            this.threePrimeBytes = StringUtil.stringToBytes(threePrimeReadOrder);
            this.fivePrimeReadOrder = fivePrimeReadOrder;
            this.fivePrimeReadOrderBytes = StringUtil.stringToBytes(fivePrimeReadOrder);
            this.fivePrime = SequenceUtil.reverseComplement(fivePrimeReadOrder);
            this.fivePrimeBytes = StringUtil.stringToBytes(this.fivePrime);
        }

        public String get3PrimeAdapter(){ return threePrime; }
        public String get5PrimeAdapter(){ return fivePrime; }
        public String get3PrimeAdapterInReadOrder(){ return threePrime; }
        public String get5PrimeAdapterInReadOrder() { return fivePrimeReadOrder; }
        public byte[] get3PrimeAdapterBytes() { return threePrimeBytes; }
        public byte[] get5PrimeAdapterBytes() { return fivePrimeBytes; }
        public byte[] get3PrimeAdapterBytesInReadOrder() { return threePrimeBytes; }
        public byte[] get5PrimeAdapterBytesInReadOrder()  { return fivePrimeReadOrderBytes; }

        public String getName() { return this.name; }

        public void setName(final String name) {
            this.name = name;
        }

        // WARNING: These methods ignore the name member!
        @Override
        public boolean equals(final Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            final TruncatedAdapterPair that = (TruncatedAdapterPair) o;

            if (!fivePrime.equals(that.fivePrime)) return false;
            if (!threePrime.equals(that.threePrime)) return false;

            return true;
        }

        @Override
        public int hashCode() {
            int result = fivePrime.hashCode();
            result = 31 * result + threePrime.hashCode();
            return result;
        }

        @Override
        public String toString() {
            return "TruncatedAdapterPair{" +
                    "fivePrimeReadOrder='" + fivePrimeReadOrder + '\'' +
                    ", threePrime='" + threePrime + '\'' +
                    ", name='" + name + '\'' +
                    '}';
        }
    }
}
