package org.broadinstitute.hellbender.utils;


import org.broadinstitute.hellbender.exceptions.GATKException;

import java.util.stream.Collectors;

/**
 * Class used for storing a list of doubles as a run length encoded histogram that compresses the data into bins spaced
 * at defined intervals.
 */
public final class Histogram {

    private final String PRINT_DELIM = ",";

    private static final double BIN_EPSILON = 0.01;


    private final double binSize;
    private final String precisionFormat;

    private MultiSet<Long> dataList = new MultiSet<>();

    /**
     * Create an empty histogram object with a default bin size of 0.1
     */
    public Histogram() {
        this.binSize = 0.1;
        this.precisionFormat = "%.1f";
    }

    /**
     * Create an empty histogram object with a with a specified bin size
     * @param binSize size of the bins to compress the data into
     */
    public Histogram(final double binSize) {
        this.binSize = binSize;
        precisionFormat = "%." + Math.round(-Math.log10(binSize)) + "f";
    }

    /**
     * Add a value to be stored in the histogram
     * @param d Data to be added to the histogram
     */
    public void add(final Double d) {
        if (d.isNaN()) {
            return;
        }
        long binKey = getBinnedValue(d);
        if (isValidBinKey(binKey)) {
            dataList.add(binKey);
        } else {
            throw new GATKException("Histogram values are suspiciously extreme.  Failed to add " + d + " to the Histogram.");
        }
    }

    /**
     * Add multiple copies of the same value into the histogram to be stored in the histogram
     * @param d Data to be added to the histogram
     * @param count number of repetitions of the value to add to the histogram
     */
    public void add(final Double d, final int count) {
        if (count < 1) {
            throw new GATKException("Cannot add non-positive counts to Histogram.");
        }
        long binKey = getBinnedValue(d);
        if (isValidBinKey(binKey)) {
            dataList.add(binKey, count);
        } else {
            throw new GATKException("Histogram values are suspiciously extreme.  Failed to add " + d + " to the Histogram.");
        }
    }

    /**
     * Combine two histogram objects that have the same bin by adding the values to each bin
     * @param h histogram to add
     */
    public void add(final Histogram h) {
        if (this.binSize != h.binSize) {
            throw new GATKException("Histogram bin sizes are mismatched -- cannot add bin size " + this.binSize + " to " + h.binSize);
        }
        this.dataList.addAll(h.dataList);
    }

    /**
     * Return the count of items in the bin corresponding to the provided value
     * @param d Value to test
     * @return Number of items in the bin mapped by the provided key
     */
    public Integer get(final double d) {
        long binKey = getBinnedValue(d);
        if (isValidBinKey(binKey)) {
            return dataList.multificity(binKey);
        } else {
            throw new GATKException("Requested value is suspiciously extreme.  Failed to retrieve " + d + " from the Histogram.");
        }
    }

    /**
     *
     * @return may be null if Histogram is empty
     */

    public Double median() {
        final long numItems = dataList.longSize();
        boolean oddNumberValues = true;
        if(numItems % 2 == 0) {
            oddNumberValues = false;
        }
        final long medianIndex = (numItems+1)/2;
        final Long[] sortedBinKeys = dataList.distinct().stream().sorted().toArray(Long[]::new);

        long leftSize = 0;
        int firstMedian = -1;
        for (int i = 0; i < sortedBinKeys.length; i++) {
            final long key = sortedBinKeys[i];
            leftSize += dataList.multificity(sortedBinKeys[i]);
            if (leftSize > medianIndex) {
                if (firstMedian < 0) {
                    return key * binSize;
                } else {
                    return binSize * (sortedBinKeys[firstMedian] + key) / 2.0;
                }
            } else if (leftSize == medianIndex) {
                if (oddNumberValues) {
                    return key * binSize;
                } else {
                    firstMedian = i;
                }
            }
        }
        return null;
    }

    private long getBinnedValue(double d) {
        return Math.round(Math.floor((d+BIN_EPSILON*binSize)/binSize)); //add a little epsilon before division so values exactly on bin boundaries will stay in the same bin
    }

    private boolean isValidBinKey(long binnedValue) {
        return binnedValue <= Integer.MAX_VALUE && binnedValue >= Integer.MIN_VALUE;
    }

    @Override
    public String toString(){
        if (dataList.isEmpty()) {
            return "" + Double.NaN;
        } else {
            return dataList.distinct().stream().sorted()
                    .map(binKey -> "" + String.format(precisionFormat, (binKey * binSize)) + PRINT_DELIM + dataList.multificity(binKey))
                    .collect(Collectors.joining(PRINT_DELIM));

        }
    }

    /**
     *
     * @return true if there is no data in this histogram
     */
    public boolean isEmpty() {
        return dataList.isEmpty();
    }
}

