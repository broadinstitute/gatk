package org.broadinstitute.hellbender.utils;


import org.broadinstitute.hellbender.exceptions.GATKException;

import java.util.Arrays;

/**
 * Class used for storing a list of doubles as a run length encoded histogram that compresses the data into bins spaced
 * at defined intervals.
 */
public class Histogram {
    private Double binSize;
    private String precisionFormat;
    private String printDelim;
    final private Double BIN_EPSILON = 0.01;

    private CompressedDataList<Integer> dataList = new CompressedDataList<>();

    /**
     * Create an empty histogram object with a default bin size of 0.1
     */
    public Histogram() {
        this.binSize = 0.1;
        precisionFormat = "%.1f";
    }

    /**
     * Create an empty histogram object with a with a specified bin size
     * @param binSize size of the bins to compress the data into
     */
    public Histogram(final Double binSize) {
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
            dataList.add((int) binKey);
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
            dataList.add((int) binKey, count);
        } else {
            throw new GATKException("Histogram values are suspiciously extreme.  Failed to add " + d + " to the Histogram.");
        }
    }

    /**
     * Combine two histogram objects that have the same bin by adding the values to each bin
     * @param h histogram to add
     */
    public void add(final Histogram h) {
        if (!this.binSize.equals(h.binSize)) {
            throw new GATKException("Histogram bin sizes are mismatched -- cannot add bin size " + this.binSize + " to " + h.binSize);
        }
        this.dataList.add(h.dataList);
    }

    /**
     * Return the count of items in the bin corresponding to the provided value
     * @param d Value to test
     * @return Number of items in the bin mapped by the provided key
     */
    public Integer get(final Double d) {
        long binKey = getBinnedValue(d);
        if (isValidBinKey(binKey)) {
            return dataList.getValueCounts().get((int) binKey);
        } else {
            throw new GATKException("Requested value is suspiciously extreme.  Failed to retrieve " + d + " from the Histogram.");
        }
    }

    /**
     *
     * @return may be null if Histogram is empty
     */

    public Double median() {
        int numItems = 0;
        for(final int count : dataList.valueCounts.values()) {
            numItems += count;
        }
        boolean oddNumberValues = true;
        if(numItems % 2 == 0) {
            oddNumberValues = false;
        }
        int medianIndex = (numItems+1)/2;

        int counter = 0;
        Double firstMedian = null;
        for(final Integer key : dataList.valueCounts.keySet()) {
            counter += dataList.valueCounts.get(key);
            if( counter > medianIndex) {
                if (firstMedian == null) {
                    return key * binSize;
                }
                else {
                    return (firstMedian+key)/2.0*binSize;
                }
            }
            if( counter == medianIndex) {
                if (oddNumberValues) {
                    return key * binSize;
                }
                else {
                    firstMedian = (double) key;
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
        printDelim = ",";
        String str = "";
        Object[] keys = dataList.valueCounts.keySet().toArray();
        if (keys.length == 0) {
            return Double.toString(Double.NaN);
        }
        Arrays.sort(keys);
        for (Object i: keys){
            if(!str.isEmpty()) {
                str += printDelim;
            }
            str+=(String.format(precisionFormat,(double)(int)i*binSize)+printDelim+dataList.valueCounts.get(i));  //key i needs to be output with specified precision
        }
        return str;
    }

    /**
     *
     * @return true if there is no data in this histogram
     */
    public boolean isEmpty() {
        return dataList.isEmpty();
    }
}

