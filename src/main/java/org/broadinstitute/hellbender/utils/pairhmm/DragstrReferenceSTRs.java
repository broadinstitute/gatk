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
            if (sequence.length  < period << 1) {
                break;
            }
            int position, remainingToMatchUnit, carryBack, rightMargin, positionPlusPeriod, resultArrayOffset;

            // We first find the right-margin enclosing the last repeat that would match
            // a prospective repeat unit within [start,end):
            rightMargin = Math.min(end + period, sequence.length) - 1;
            final int rightMostStart = position = rightMargin - period;
            remainingToMatchUnit = period;
            carryBack = 1;
            for (; rightMargin < sequence.length; rightMargin++) {
                        if (sequence[position++] != sequence[rightMargin]) {
                    break;
                } else if (--remainingToMatchUnit == 0) {
                    carryBack++;
                    remainingToMatchUnit = period;
                }
            }

            // then we work our way backwards carrying on number of conseqcutive matching units
            // and updating.
            // carryBack and matchesRun has been updated in the forward pass so they correspond to the value
            // at rightMostStart.
            if ((resultArrayOffset = rightMostStart - start) >= 0 && carryBack > repeats[resultArrayOffset]) {
                repeats[resultArrayOffset] = carryBack;
                periods[resultArrayOffset] = period;
            }

            boolean inZone = false;
            for (position = rightMostStart - 1, positionPlusPeriod = position + period,
                 resultArrayOffset = position - start; position >= start; position--, positionPlusPeriod--, resultArrayOffset--) {
                if (sequence[position] == sequence[positionPlusPeriod]) {
                    if (--remainingToMatchUnit == 0) { // have we matched yet another unit of length period?
                        ++carryBack;
                        remainingToMatchUnit = period;
                    }
                    if ((inZone |= position < end) && carryBack > repeats[resultArrayOffset]) {
                        repeats[resultArrayOffset] = carryBack;
                        periods[resultArrayOffset] = period;
                    }
                } else {
                    carryBack = 1;
                    remainingToMatchUnit = period;
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
