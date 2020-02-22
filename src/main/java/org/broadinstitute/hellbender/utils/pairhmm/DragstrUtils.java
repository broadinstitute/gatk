package org.broadinstitute.hellbender.utils.pairhmm;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.aeonbits.owner.util.Collections;
import org.broadinstitute.hellbender.utils.MathUtils;

import java.util.Arrays;
import java.util.Collection;

public class DragstrUtils {

    public static STRSequenceAnalyzer repeatPeriodAndCounts(final int maxSequenceLength, final int maxPeriod) {
        return new STRSequenceAnalyzer(maxSequenceLength, maxPeriod);
    }

    public static STRSequenceAnalyzer repeatPeriodAndCounts(final byte[] sequence, final int maxPeriod) {
        final STRSequenceAnalyzer result = new STRSequenceAnalyzer(sequence.length, maxPeriod);
        result.load(sequence);
        return result;
    }

    public static STRSequenceAnalyzer repeatPeriodAndCounts(final byte[] sequence, final int start, final int stop, final int maxPeriod) {
        final STRSequenceAnalyzer result = new STRSequenceAnalyzer(sequence.length, maxPeriod);
        result.load(sequence, start, stop);
        return result;
    }

    public static Collection<? extends VCFHeaderLine> vcfHeaderLines() {
        return Collections.list(new VCFInfoHeaderLine(DragstrConstants.DRAGSTRINFO_KEY, 2, VCFHeaderLineType.Integer, "Indicates the period and repeat count"),
                                new VCFInfoHeaderLine(DragstrConstants.DRAGSTRPARAMS_KEY, 3, VCFHeaderLineType.Float, "Parameeters used (GOP, GCP, API)"));

    }

    public static VariantContext annotate(VariantContext annotatedCall, final DragstrParams dragstrParams, STRSequenceAnalyzer dragstrs, int offset, int ploidy, double snpHeterozygosity) {
        final VariantContextBuilder builder = new VariantContextBuilder(annotatedCall);
        final int period = dragstrs.mostRepeatedPeriod(offset);
        final int repeats = dragstrs.numberOfMostRepeats(offset);
        final double gop = dragstrParams.gop(period, repeats);
        final double gcp = dragstrParams.gcp(period, repeats);
        final double api = dragstrParams.api(period, repeats);
        builder.attribute(DragstrConstants.DRAGSTRINFO_KEY, new int[] {period, repeats});
        builder.attribute(DragstrConstants.DRAGSTRPARAMS_KEY, new String[] {String.format("%.1f", gop),String.format("%.1f", gcp), String.format("%.1f", api)});
        return builder.make();
    }

    public static class STRSequenceAnalyzer {
        private int start;
        private int end;
        private final int[][] repeatsByPeriodAndPosition;
        private final int[] periodWithMostRepeats;
        private final int maxPeriod;

        private STRSequenceAnalyzer(final int maxSequenceLength, final int maxPeriod) {
            repeatsByPeriodAndPosition = new int[maxPeriod][maxSequenceLength];
            this.maxPeriod = maxPeriod;
            this.periodWithMostRepeats = new int[maxSequenceLength];
        }

        public int numberOfRepeats(final int position, final int period) {
            if (period <= 0 || period > maxPeriod) {
                return 0;
            } else if (position < start || position >= end) {
                throw new IllegalArgumentException("cannot query outside requested boundaries");
            } else {
                return repeatsByPeriodAndPosition[period - 1][position];
            }
        }

        public int mostRepeatedPeriod(final int position) {
            if (position >= start && position < end) {
                return periodWithMostRepeats[position];
            } else {
                throw new IllegalArgumentException("cannot query outside requested boundaries");
            }
        }

        public int numberOfMostRepeats(final int position) {
            if (position >= start && position < end) {
                return repeatsByPeriodAndPosition[periodWithMostRepeats[position] - 1][position];
            } else {
                throw new IllegalArgumentException("cannot query outside requested boundaries");
            }
        }

        public void load(final byte[] sequence) {
            load(sequence, 0, sequence.length);
        }

        public void load(final byte[] sequence, final int start, final int end) {
            if (sequence.length > repeatsByPeriodAndPosition[0].length) {
                throw new IllegalArgumentException("input sequence is too long");
            } else if (end < start || start < 0 || end > sequence.length) {
                throw new IndexOutOfBoundsException("bad indexes " + start  + " " + end + " " + sequence.length);
            } else if (start >= end) {
                this.start = start; this.end = end;
                return;
            }

            this.start = start;
            this.end = end;
            // periodIndex == periodLength - 1 since in Java indexes start with 0.

            loadPeriodOne(sequence, start, end); // simpler and faster code for period-length == 0.
            for (int periodIndex = 1; periodIndex < maxPeriod; periodIndex++) {
                final int periodLength = periodIndex + 1;
                final int[] runLength = repeatsByPeriodAndPosition[periodIndex];
                if (sequence.length < periodLength) {
                    Arrays.fill(runLength, 0, sequence.length, 0);
                    continue;
                }
                int position, matchedCycles, prevValue, cycleIndex, leftMargin, rightMargin;

                // we calculate the margins, first right-margin, by looking for
                // the first mismatch after then sequence of interest.
                int stop;
                if (end < sequence.length - periodLength) {
                    stop = sequence.length - periodLength;
                    for (position = end - 1, rightMargin = Math.min(position + periodLength, sequence.length); position < stop; position++, rightMargin++) {
                        if (sequence[position] != sequence[rightMargin]) {
                            break;
                        }
                    }
                } else {
                    rightMargin = sequence.length;
                }

                // Calculate backward the repeat run lengths from a position forward.
                // so that runLength[i] would contain the number of repetitions of seq[i .. i + period) in seq[i .. end_of_seq]
                // First We set period - 1 last positions to 0 as the units will fall off the end of read.
                for (position = rightMargin - 1, cycleIndex = 1; cycleIndex < periodLength; position--, cycleIndex++) {
                    runLength[position] = 0;
                }

                // Then the main loop.
                // We basically go backward checking what was the number of reported repeats downstream
                // and increasing as we get to a number exact matched cycles that equal to the period (length).
                // When we found a mismatch we reset that count to 0 and the run length to 1 (the minimum run-length).
                prevValue = runLength[position--] = 1; //prevValue holds the num of repeats reported in the previous (+1) position.
                stop = start - periodLength;
                for (leftMargin = position, position += periodLength, matchedCycles = 0; leftMargin >= 0; leftMargin--, position--) {
                    if (sequence[leftMargin] == sequence[position]) { // we keep matching repeat unit bases.
                        if (++matchedCycles == periodLength) { // we go a new full repeat matched so the run length increases:
                            prevValue = runLength[leftMargin] = prevValue + 1;
                            matchedCycles = 0; // we reset the match-run-length to 0.
                        } else { // we simply copy the run length from the +1 position base.
                            runLength[leftMargin] = prevValue;
                        }
                    } else if (leftMargin >= stop){ // we bump into a mismatch that end the run.
                        prevValue = runLength[leftMargin] = 1;
                        matchedCycles = 0; // we reset the match-run-length to 0.
                    } else {
                        break;
                    }
                }
                leftMargin++;

                // propagate forward the total run-length over to the other repeats in the run.
                // we do it per cycle in the period so that we first deal with repeat units whose offset
                // is zero respect the beginning of the sequence, the2 1 base, then 2 etc.
                for (cycleIndex = 0; cycleIndex < periodLength; cycleIndex++) {
                  // The left most repeated unit runLength[i] contains the actual run length for all
                  // the units to the right. We copy that value forward to the other run-length units.
                  // We do this by iterating over consecutive repeat runs.
                  for (position = leftMargin + cycleIndex; position < rightMargin; position += periodLength) {
                    final int totalRunLength = runLength[position];
                    for (int repeatInRun = 1; repeatInRun < totalRunLength; repeatInRun++) {
                        runLength[position += periodLength] = totalRunLength;
                    }
                  }
                }

                // Now we calculate the max repeat length that overlaps any given position.

                // we skip period == 1 (periodIndex == 0) since is already resolved.

                // for period == 2 the code can be simplified a bit since it only requires the combination of the
                // current value and the previous one
                if (periodLength == 2 && end > 2) { // end > 2 to avoid a leftIndex = -1 condition.
                    final int[] twoPeriodValues = repeatsByPeriodAndPosition[1]; // this way we avoid repeated indirection two the period 2 array.
                    int rightIndex = end - 1,
                        leftIndex = rightIndex - 1,
                        rightValue = twoPeriodValues[rightIndex],
                        leftValue;
                    do {
                        leftValue = twoPeriodValues[leftIndex];
                        if (leftValue > rightValue) {
                            twoPeriodValues[rightIndex] = leftValue;
                        }
                        rightValue = leftValue;
                        rightIndex--;
                        leftIndex--;
                    } while (leftIndex >= leftMargin);
                } else if (periodLength > 2 && sequence.length > 2) {
                    // for period 3 or above we could use a special heap to get in k log k the longest length in a window
                    // but since in practice the max perior is something like 8 or 10 bases it seems a bit of an overkill.
                        final int[] periodValues = repeatsByPeriodAndPosition[periodIndex];
                        // we fill the expected trailing 0 with the count for the first non-zero which is exactly
                        // at length - periodLength (e.q. length - periodIndex - 1).
                        int windowEnd, windowStart, maxInWindow;
                        windowStart = rightMargin - periodLength;
                        windowEnd = rightMargin;
                        maxInWindow = periodValues[windowStart]; // guaranteed to be the the best as other are always 0.
                        while (windowEnd > 1) { // not > 0 since the very first value does not need to be changed.
                            int valueOut = periodValues[--windowEnd];
                            int valueIn = --windowStart < leftMargin ? -1 : periodValues[windowStart];
                            periodValues[windowEnd] = maxInWindow;
                            if (valueIn < valueOut && valueOut >= maxInWindow) {
                                // here we pay a bit of a penalty to make sure that we get the maxInWindow value correct since
                                // we are removing a value that is the same (or larger) than the current maxInWindow ant the value
                                // coming in is not larger than that.
                                maxInWindow = MathUtils.arrayMax(periodValues, windowStart < leftMargin ? leftMargin : windowStart, windowEnd, -1);
                            } else if (valueIn > maxInWindow) {
                                // the value coming in is larger than the current max so the maxInWindow update is trivial.
                                maxInWindow = valueIn;
                            }
                        }
                }
            }

            // finally we update the periodWithMostRepeats array:

            Arrays.fill(periodWithMostRepeats, start, end, 1);
            final int[] mostRepeats = repeatsByPeriodAndPosition[0].clone();
            for (int periodIndex = 1; periodIndex < maxPeriod; periodIndex++) {
                final int periodLength = periodIndex + 1;
                final int[] periodValues = repeatsByPeriodAndPosition[periodIndex];
                for (int position = start; position < end; position++) {
                    final int repeats = periodValues[position];
                    if (repeats > mostRepeats[position]) {
                        mostRepeats[position] = repeats;
                        periodWithMostRepeats[position] = periodLength;
                    }
                }
            }
        }

        private void loadPeriodOne(final byte[] sequence, final int start, final int end) {
            final int[] runLengths = repeatsByPeriodAndPosition[0];
            byte last = sequence[end - 1];

            // backward phase:
            int rightMargin, position;
            for (rightMargin = end; rightMargin < sequence.length && sequence[rightMargin] == last; rightMargin++);
            int carryBack = rightMargin - end;
            for (position = end - 1; position >= start; position--) {
                final byte next = sequence[position];
                runLengths[position] = next == last ? ++carryBack : (carryBack = 1);
                last = next;
            }
            // forward phase:
            int leftMargin;
            last = sequence[start];
            for (leftMargin = start - 1; leftMargin >= 0 && sequence[leftMargin] == last; leftMargin--);
            int carryForward = start - leftMargin - 1;
            for (position = start; position < end; position++) {
                final byte next = sequence[position];
                if (next == last) {
                    runLengths[position] += carryForward++;
                } else {
                    carryForward = 1;
                }
                last = next;
            }
        }
    }

}
