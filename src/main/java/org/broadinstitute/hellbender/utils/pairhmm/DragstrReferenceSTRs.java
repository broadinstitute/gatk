package org.broadinstitute.hellbender.utils.pairhmm;

import org.apache.commons.lang.ArrayUtils;

import java.util.Arrays;

public class DragstrReferenceSTRs {
    private final int start;
    private final int end;
    private final byte[] bases;
    private final int[] period;
    private final int[] forwardRepeats;

    private DragstrReferenceSTRs(final byte[] bases, final int start, final int end, final int[] period, final int[] repeats) {
        this.start = start;
        this.end = end;
        this.bases = bases;
        this.period = period;
        this.forwardRepeats = repeats;
    }

    public int repeatLength(final int position) {
        int result = lookup(position, forwardRepeats);
        final int period = this.period[position - start];
        for (int i = position - 1, j = position + period, k = period; i >= 0; i--) {
            if (bases[i] != bases[--j]) {
                break;
            } else if (--k == 0) { // here we count complete unit matches (k reaches 0).
                k = period;
                result++;
            }
        }
        return result;
    }

    public byte[] repeatUnit(final int position) {
        final int length = lookup(position, period);
        return Arrays.copyOfRange(bases, position, position + length);
    }

    public String repeatUnitAsString(final int position) {
        return new String(repeatUnit(position));
    }

    public int period(final int position) {
        return lookup(position, period);
    }

    private int lookup(final int position, final int[] array) {
        final int offset = position - start;
        if (offset >= 0 && position <= end) {
            return array[offset];
        } else {
            throw new IllegalArgumentException("postion " + position + " is outside bounds");
        }
    }


    public static DragstrReferenceSTRs of(final byte[] sequence, final int start, final int end, final int maxPeriod) {
        if (end < start || start < 0 || end > sequence.length) {
            throw new IndexOutOfBoundsException("bad indexes " + start  + " " + end + " " + sequence.length);
        } else if (start >= end) {
            return new DragstrReferenceSTRs(sequence, start, end, ArrayUtils.EMPTY_INT_ARRAY, ArrayUtils.EMPTY_INT_ARRAY);
        }
        final int[] repeats = processPeriodOne(sequence, start, end);
        final int[] periods = new int[end - start];
        Arrays.fill(periods, 1);
        for (int period = 2; period <= maxPeriod; period++) {
            if (sequence.length < period << 1) {
                break;
            }
            int position, matchesRun, carryBack, rightMargin;

            // We first find the right-margin enclosing the last repeat that would match
            // a prospective repeat unit within [start,end):
            rightMargin = Math.min(end + period, sequence.length) - 1;
            position = rightMargin - period;
            for (; rightMargin < sequence.length; rightMargin++) {
                if (sequence[position] != sequence[rightMargin]) {
                    break;
                }
            }

            // then we work our way backwards carrying on number of conseqcutive matching units
            // and updating.

            // Is a bit silly that we revisit each position between margin and end since the repeats&periods array
            // would be updated but it simplifies the initialization of these three vars just below. Could be optimized.
            position = rightMargin - period - 1;
            carryBack = 1;
            matchesRun = 0;
            boolean inZone = false;
            for (int compareWithPosition = position + period, offset = position - start; position >= start; position--, compareWithPosition--, offset--) {
                if (sequence[position] == sequence[compareWithPosition]) {
                    if (++matchesRun == period) { // have we matched yet another unit of length period?
                        ++carryBack;
                        matchesRun = 0;
                    }
                    if ((inZone |= position < end) && carryBack > repeats[offset]) {
                        repeats[offset] = carryBack;
                        periods[offset] = period;
                    }
                } else {
                    carryBack = 1;
                    matchesRun = 0;
                }
            }
        }
        return new DragstrReferenceSTRs(sequence, start, end, periods, repeats);
    }

    private static int[] processPeriodOne(final byte[] sequence, final int start, final int end) {
        final int[] repeats = new int[end - start];
        byte last = sequence[end - 1];

        // backward phase:
        int rightMargin, position;
        for (rightMargin = end; rightMargin < sequence.length && sequence[rightMargin] == last; rightMargin++);
        int carryBack = rightMargin - end;
        int offset;
        for (position = end - 1, offset = position - start; position >= start; position--, offset--) {
            final byte next = sequence[position];
            repeats[offset] = next == last ? ++carryBack : (carryBack = 1);
            last = next;
        }
        return repeats;
    }
}
