package org.broadinstitute.hellbender.utils.pairhmm;

import org.apache.commons.lang.ArrayUtils;

import java.util.Arrays;

public class DragstrReadSTRAnalyzer {
    private final int[][] repeatsByPeriodAndPosition;
    private final int[] periodWithMostRepeats;
    private final int maxPeriod;
    /**
     * Length of the base sequence
     */
    private int seqLength;

    // array use as a temporary variable to hold a small heap.
    private transient int[] heap;
    // array use as a temporary variable to hold run-lengths.
    private transient int[] tempRunLengths;

    DragstrReadSTRAnalyzer(final int maxSequenceLength, final int maxPeriod) {
        repeatsByPeriodAndPosition = new int[maxPeriod][maxSequenceLength];
        this.maxPeriod = maxPeriod;
        this.periodWithMostRepeats = new int[maxSequenceLength];
        this.heap = new int[maxPeriod + 1];
        this.tempRunLengths = new int[maxSequenceLength];
    }

    public int numberOfRepeats(final int position, final int period) {
        if (period <= 0 || period > maxPeriod) {
            return 0;
        } else if (position < 0 || position >= seqLength) {
            throw new IllegalArgumentException("cannot query outside requested boundaries");
        } else {
            return repeatsByPeriodAndPosition[period - 1][position];
        }
    }

    public int mostRepeatedPeriod(final int position) {
        if (position >= 0 && position < seqLength) {
            return periodWithMostRepeats[position];
        } else {
            throw new IllegalArgumentException(String.format("cannot query outside requested boundaries [%d, %d): %d", 0, seqLength, position));
        }
    }

    public int numberOfMostRepeats(final int position) {
        if (position >= 0 && position < seqLength) {
            return repeatsByPeriodAndPosition[periodWithMostRepeats[position] - 1][position];
        } else {
            throw new IllegalArgumentException("cannot query outside requested boundaries");
        }
    }

    public void load(final byte[] sequence) {
        if (sequence.length > repeatsByPeriodAndPosition[0].length) {
            throw new IllegalArgumentException("input sequence is too long");
        }
        seqLength = sequence.length;
        // periodIndex == periodLength - 1 since in Java indexes start with 0.

        calculateRepeatsForPeriodOne(sequence); // simpler and faster code for period-length == 0.
        for (int period = 2; period <= maxPeriod; period++) {
            calculateRepeatsForPeriodTwoAndAbove(period, sequence);
        }

        // finally we update the periodWithMostRepeats array
        // we the best period at each position.
        Arrays.fill(periodWithMostRepeats, 0, seqLength, 1);
        final int[] mostRepeats = repeatsByPeriodAndPosition[0].clone();
        for (int periodIndex = 1; periodIndex < maxPeriod; periodIndex++) {
            final int periodLength = periodIndex + 1;
            final int[] periodValues = repeatsByPeriodAndPosition[periodIndex];
            for (int position = 0; position < seqLength; position++) {
                final int repeats = periodValues[position];
                if (repeats > mostRepeats[position]) {
                    mostRepeats[position] = repeats;
                    periodWithMostRepeats[position] = periodLength;
                }
            }
        }
    }

    private void calculateRepeatsForPeriodTwoAndAbove(final int period, final byte[] sequence) {
        // very small sequence (less than period) then repeats are all zero:2
        if (seqLength < period) {
            Arrays.fill(tempRunLengths, 0, seqLength, 0);
            return;
        }

        int position, matchedCycles, cycleIndex, positionPlusPeriod;

        // suffixes of zero run-lengths at the end of the sequence (period - 1).
        for (position = seqLength - 1, cycleIndex = period; cycleIndex > 1; position--, cycleIndex--) {
            tempRunLengths[position] = 0;
        }

        // backward phase.
        // at the end of this phase runLength[x] will have the total number of equal repeats downstream with length starting at position X
        // inclusive.
        int carryBack = tempRunLengths[position--] = 1; //prevValue holds the num of repeats reported in the previous (+1) position.
        for (positionPlusPeriod = position + period, matchedCycles = 0; position >= 0; position--, positionPlusPeriod--) {
            if (sequence[position] == sequence[positionPlusPeriod]) { // we keep matching repeat unit bases.
                if (++matchedCycles == period) { // we go a new full repeat matched so the run length increases:
                    tempRunLengths[position] = ++carryBack;
                    matchedCycles = 0; // we reset the match-run-length to 0.
                } else { // we simply copy the run length from the +1 position base.
                    tempRunLengths[position] = carryBack;
                }
            } else { // we bump into a mismatch that end the run.
                carryBack = tempRunLengths[position] = 1;
                matchedCycles = 0; // we reset the match-run-length to 0.
            }
        }
        // Now we propagate forward the number of equivalent repeats up stream:
        // So at the end runLength[X] == the repeats (up and down-stream) of the
        // unit of length period that starts at X.
        // The outside for loop iterates along different repeat unit start offsets (here named cycles):
        for (cycleIndex = 0; cycleIndex < period; cycleIndex++) {
            // The left most repeated unit runLength[i] contains the actual run length for all
            // the units to the right. We copy that value forward to the other run-length units.
            // We do this by iterating over consecutive repeat runs.
            for (position = cycleIndex; position < seqLength; position += period) {
                final int totalRunLength = tempRunLengths[position];
                for (int repeatInRun = 1; repeatInRun < totalRunLength; repeatInRun++) {
                    tempRunLengths[position += period] = totalRunLength;
                }
            }
        }

        // NOTE: this line below is added solely to replicate an oddity of DRAGEN's algorithm that discounts one repeat
        // only at the beginning of the sequence for period 2 or above (e.g. ^CACATG would yield periods 122211 repeats
        // 12221 when the natural/intuitive solution is period 222211 reps 222211). This is true for the original Matlab
        // and latest DRAGEN so we must add this here:
        tempRunLengths[0]--; // No need to bother to check than repeats decreases to less than 1 since if so, period == 1 will win anyways.

        // Finally we need to propagate the best run-lengths to neighboor positions
        // so that a position has the maximum run-length of all the possible
        // units that contain the position + the unit that follow up-stream.
        final int[] finalRunLengths = repeatsByPeriodAndPosition[period - 1];
        final int heapSize = period + 1;
        Arrays.fill(heap, 0, heapSize, 0);
        int currentMax = heap[0] = tempRunLengths[0];
        // the first period length prefix is trivial:
        // a position's max run-lenght is the largest seen so far and the following position's
        // We use this opportunity to populate the
        // run-length max-heap
        for (position = 0; position < period; ) {
            finalRunLengths[position] = currentMax = Math.max(currentMax, heap[++position] = tempRunLengths[position]);
            fixHeap(position, heapSize);
        }
        for (int outPosition = 0; ++position < seqLength; ) {
            final int valueOut = tempRunLengths[outPosition++];
            final int valueIn = tempRunLengths[position];
            if (valueIn != valueOut) {
                final int fixHeapIdx = ArrayUtils.indexOf(heap, valueOut); // O(n) search although usually n is very small.
                heap[fixHeapIdx] = valueIn;
                fixHeap(fixHeapIdx, heapSize);
                currentMax = heap[0];
            }
            finalRunLengths[position - 1] = currentMax;
        }
        finalRunLengths[position - 1] = currentMax; // we need to do the last one.
    }

    private void fixHeap(final int idx, final int heapSize) {
        if (idx == 0 || heap[(idx - 1) >> 1] > heap[idx]) {
            fixHeapDown(idx, heapSize);
        } else {
            fixHeapUp(idx);
        }
    }

    private void fixHeapUp(int idx) {
        final int value = heap[idx];
        do {
            final int upIdx = (idx - 1) >> 1;
            final int upValue = heap[upIdx];
            if (upValue >= value) {
                break;
            } else {
                heap[idx] = upValue;
                idx = upIdx;
            }
        } while (idx > 0);
        heap[idx] = value;
    }

    // Method useful when debugging, keep it around just in case.

    /**
     * Check that the content of the heap is correct given its size.
     *
     * Used for debugging, it might be useful in the future; better to keep it
     * around.
     */
    @SuppressWarnings("unused")
    private boolean checkHeap(final int heapSize) {
        for (int i = 1; i < heapSize; i++) {
            if (heap[(i - 1) >> 1] < heap[i]) {
                return false;
            }
        }
        return true;
    }

    private void fixHeapDown(int idx, final int heapSize) {
        final int value = heap[idx];
        while (true) {
            final int idxRight = (idx + 1) << 1;
            final int rightValue = idxRight < heapSize ? heap[idxRight] : -1;
            final int idxLeft = idxRight - 1;
            final int leftValue = idxLeft < heapSize ? heap[idxLeft] : -1;
            if (rightValue > value) {
                if (rightValue > leftValue) {
                    heap[idx] = rightValue;
                    idx = idxRight;
                } else {
                    heap[idx] = leftValue;
                    idx = idxLeft;
                }
            } else if (leftValue > value) {
                heap[idx] = leftValue;
                idx = idxLeft;
            } else {
                break;
            }
        }
        heap[idx] = value;
    }

    /**
     * I faster and simpler implementation of {@link #calculateRepeatsForPeriodTwoAndAbove(int, byte[])}
     * for period == 1.
     * <p>
     *     Avoids the unecessary matched cycles book-keeping and the use of  very trivial heap.
     * </p>
     * @param sequence
     */
    private void calculateRepeatsForPeriodOne(final byte[] sequence) {
        final int[] runLengths = repeatsByPeriodAndPosition[0];
        final int rightMargin = seqLength - 1;
        byte last = sequence[rightMargin];

        // backward phase:
        int carryBack = runLengths[rightMargin] = 1;
        for (int position = rightMargin - 1; position >= 0; position--) {
            final byte next = sequence[position];
            runLengths[position] = next == last ? ++carryBack : (carryBack = 1);
            last = next;
        }
        // forward phase:
        // last = sequence[0]; // already true.
        int prevRunLength = runLengths[0];
        for (int position = 1; position <= rightMargin; position++) {
            final byte next = sequence[position];
            if (next == last) {
                runLengths[position] = prevRunLength;
            } else {
                final int thisRunLength = runLengths[position];
                if (prevRunLength < thisRunLength) {
                    runLengths[position - 1] = thisRunLength;
                }
                last = next;
                prevRunLength = thisRunLength;
            }
        }
    }
}
