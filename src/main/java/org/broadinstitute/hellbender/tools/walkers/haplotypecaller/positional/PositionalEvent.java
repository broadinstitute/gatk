package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.positional;

import org.broadinstitute.gatk.nativebindings.smithwaterman.SWParameters;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;

/**
 * A primitive event (a substitution of one base, an insertion or a deletion) in the coordinates of a positional
 * assembly's reference window: {@code ref} is replaced by {@code alt} at window offset {@code position}.  An insertion
 * has an empty {@code ref} and goes before {@code position}; a deletion has an empty {@code alt}.
 */
record PositionalEvent(int position, String ref, String alt) implements Comparable<PositionalEvent> {

    /** offset of the first reference base after the event */
    int end() {
        return position + ref.length();
    }

    boolean isSubstitution() {
        return !ref.isEmpty() && ref.length() == alt.length();
    }

    /**
     * This event left-aligned as far as the reference allows: the same for every position an indel can take in a
     * tandem repeat, whatever other events stopped a read's copy of it short of the leftmost one.
     */
    PositionalEvent leftAligned(final String reference) {
        if (ref.isEmpty()) {
            return leftAlignedInsertion(reference, position, alt, 0);
        } else if (alt.isEmpty()) {
            return leftAlignedDeletion(reference, position, ref.length(), 0);
        }
        return this;
    }

    @Override
    public int compareTo(final PositionalEvent other) {
        final int byPosition = Integer.compare(position, other.position);
        if (byPosition != 0) {
            return byPosition;
        }
        final int byRef = ref.compareTo(other.ref);
        return byRef != 0 ? byRef : alt.compareTo(other.alt);
    }

    /**
     * Apply non-overlapping events to {@code sequence}, which starts at window offset {@code offset}.
     * @return the edited sequence, or {@code null} if events overlap or fall outside the sequence
     */
    static String apply(final String sequence, final int offset, final Collection<PositionalEvent> events) {
        final List<PositionalEvent> sorted = new ArrayList<>(events);
        Collections.sort(sorted);
        final StringBuilder result = new StringBuilder(sequence.length() + 16);
        int next = 0;
        for (final PositionalEvent event : sorted) {
            final int start = event.position - offset;
            if (start < next || start + event.ref.length() > sequence.length()) {
                return null;
            }
            result.append(sequence, next, start).append(event.alt);
            next = start + event.ref.length();
        }
        return result.append(sequence, next, sequence.length()).toString();
    }

    /**
     * Turn the replacement of {@code reference[position, position + ref.length())} (where {@code ref} must equal those
     * reference bases) by {@code alt} into primitive events: trim common bases, decompose with an affine-gap alignment
     * scored as haplotypes are aligned to the reference, and left-align indels, never past {@code floor} or the
     * preceding event, so that the events never overlap.
     */
    static List<PositionalEvent> normalize(final String reference, final int position, final String ref, final String alt,
                                           final int floor, final SWParameters scores) {
        int refEnd = ref.length(), altEnd = alt.length();
        while (refEnd > 0 && altEnd > 0 && ref.charAt(refEnd - 1) == alt.charAt(altEnd - 1)) {
            refEnd--;
            altEnd--;
        }
        int prefix = 0;
        while (prefix < refEnd && prefix < altEnd && ref.charAt(prefix) == alt.charAt(prefix)) {
            prefix++;
        }
        final String trimmedRef = ref.substring(prefix, refEnd);
        final String trimmedAlt = alt.substring(prefix, altEnd);
        final int start = position + prefix;
        if (trimmedRef.isEmpty() && trimmedAlt.isEmpty()) {
            return Collections.emptyList();
        } else if (trimmedRef.length() == 1 && trimmedAlt.length() == 1) {
            return Collections.singletonList(new PositionalEvent(start, trimmedRef, trimmedAlt));
        } else if (trimmedAlt.isEmpty()) {
            return Collections.singletonList(leftAlignedDeletion(reference, start, trimmedRef.length(), floor));
        } else if (trimmedRef.isEmpty()) {
            return Collections.singletonList(leftAlignedInsertion(reference, start, trimmedAlt, floor));
        }

        final List<PositionalEvent> events = new ArrayList<>();
        int eventFloor = floor;
        for (final AlignmentOperation op : align(trimmedRef, trimmedAlt, scores)) {
            final PositionalEvent event;
            if (op.type == AlignmentOperation.Type.MISMATCH) {
                event = new PositionalEvent(start + op.refIndex, trimmedRef.substring(op.refIndex, op.refIndex + 1),
                        trimmedAlt.substring(op.altIndex, op.altIndex + 1));
                eventFloor = event.position + 1;
            } else if (op.type == AlignmentOperation.Type.DELETION) {
                event = leftAlignedDeletion(reference, start + op.refIndex, op.length, eventFloor);
                eventFloor = event.position + op.length;
            } else {
                event = leftAlignedInsertion(reference, start + op.refIndex,
                        trimmedAlt.substring(op.altIndex, op.altIndex + op.length), eventFloor);
                eventFloor = event.position + 1;
            }
            events.add(event);
        }
        return events;
    }

    private static PositionalEvent leftAlignedDeletion(final String reference, final int position, final int length, final int floor) {
        int start = position;
        while (start > floor && reference.charAt(start - 1) == reference.charAt(start + length - 1)) {
            start--;
        }
        return new PositionalEvent(start, reference.substring(start, start + length), "");
    }

    private static PositionalEvent leftAlignedInsertion(final String reference, final int position, final String inserted, final int floor) {
        int start = position;
        String bases = inserted;
        while (start > floor && reference.charAt(start - 1) == bases.charAt(bases.length() - 1)) {
            bases = reference.charAt(start - 1) + bases.substring(0, bases.length() - 1);
            start--;
        }
        return new PositionalEvent(start, "", bases);
    }

    private record AlignmentOperation(Type type, int refIndex, int altIndex, int length) {
        enum Type { MISMATCH, DELETION, INSERTION }
    }

    private static final int IMPOSSIBLE = Integer.MIN_VALUE / 4;
    private static final int MATCH_STATE = 0, DELETION_STATE = 1, INSERTION_STATE = 2;

    /**
     * Global affine-gap alignment of two non-empty strings.  Ties prefer a match state, then a deletion, then an
     * insertion.  Returns the mismatches and the gaps, with runs of single-base gaps merged into indels.
     */
    private static List<AlignmentOperation> align(final String ref, final String alt, final SWParameters scores) {
        final int n = ref.length(), m = alt.length(), width = m + 1;
        final int match = scores.getMatchValue(), mismatch = scores.getMismatchPenalty(),
                open = scores.getGapOpenPenalty(), extend = scores.getGapExtendPenalty();
        // score and back-pointer (the previous state) of each state of each cell
        final int[] score = new int[(n + 1) * width * 3];
        final byte[] back = new byte[(n + 1) * width * 3];
        java.util.Arrays.fill(score, IMPOSSIBLE);
        score[0] = 0;
        for (int i = 1; i <= n; i++) {
            score[3 * (i * width) + DELETION_STATE] = open + (i - 1) * extend;
            back[3 * (i * width) + DELETION_STATE] = (byte) (i > 1 ? DELETION_STATE : MATCH_STATE);
        }
        for (int j = 1; j <= m; j++) {
            score[3 * j + INSERTION_STATE] = open + (j - 1) * extend;
            back[3 * j + INSERTION_STATE] = (byte) (j > 1 ? INSERTION_STATE : MATCH_STATE);
        }
        for (int i = 1; i <= n; i++) {
            final char refBase = ref.charAt(i - 1);
            for (int j = 1; j <= m; j++) {
                final int cell = 3 * (i * width + j);
                final int diagonal = 3 * ((i - 1) * width + j - 1), up = 3 * ((i - 1) * width + j), left = 3 * (i * width + j - 1);
                final int s = refBase == alt.charAt(j - 1) ? match : mismatch;
                int best = best(score[diagonal], score[diagonal + 1], score[diagonal + 2]);
                score[cell] = score[diagonal + best] + s;
                back[cell] = (byte) best;
                best = best(score[up] + open, score[up + 1] + extend, score[up + 2] + open);
                score[cell + 1] = best == MATCH_STATE ? score[up] + open : (best == DELETION_STATE ? score[up + 1] + extend : score[up + 2] + open);
                back[cell + 1] = (byte) best;
                best = best(score[left] + open, score[left + 1] + open, score[left + 2] + extend);
                score[cell + 2] = best == MATCH_STATE ? score[left] + open : (best == DELETION_STATE ? score[left + 1] + open : score[left + 2] + extend);
                back[cell + 2] = (byte) best;
            }
        }

        final List<AlignmentOperation> reversed = new ArrayList<>();
        int i = n, j = m;
        final int last = 3 * (n * width + m);
        int state = best(score[last], score[last + 1], score[last + 2]);
        while (i > 0 || j > 0) {
            final int previous = back[3 * (i * width + j) + state];
            if (state == MATCH_STATE) {
                if (ref.charAt(i - 1) != alt.charAt(j - 1)) {
                    reversed.add(new AlignmentOperation(AlignmentOperation.Type.MISMATCH, i - 1, j - 1, 1));
                }
                i--;
                j--;
            } else if (state == DELETION_STATE) {
                reversed.add(new AlignmentOperation(AlignmentOperation.Type.DELETION, i - 1, j, 1));
                i--;
            } else {
                reversed.add(new AlignmentOperation(AlignmentOperation.Type.INSERTION, i, j - 1, 1));
                j--;
            }
            state = previous;
        }
        Collections.reverse(reversed);

        final List<AlignmentOperation> merged = new ArrayList<>();
        for (final AlignmentOperation op : reversed) {
            final AlignmentOperation previous = merged.isEmpty() ? null : merged.get(merged.size() - 1);
            final boolean extendsPrevious = previous != null && op.type != AlignmentOperation.Type.MISMATCH && previous.type == op.type &&
                    (op.type == AlignmentOperation.Type.DELETION ? previous.refIndex + previous.length == op.refIndex
                            : previous.altIndex + previous.length == op.altIndex && previous.refIndex == op.refIndex);
            if (extendsPrevious) {
                merged.set(merged.size() - 1, new AlignmentOperation(previous.type, previous.refIndex, previous.altIndex, previous.length + 1));
            } else {
                merged.add(op);
            }
        }
        return merged;
    }

    /** index of the highest of three scores, the lowest index on ties */
    private static int best(final int a, final int b, final int c) {
        int best = 0, bestScore = a;
        if (b > bestScore) {
            best = 1;
            bestScore = b;
        }
        return c > bestScore ? 2 : best;
    }
}
