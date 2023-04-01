package org.broadinstitute.hellbender.utils.read;

import org.jetbrains.annotations.NotNull;

/**
 * A CharSequence that is a view of the reverse-complement of another sequence.
 */
public final class SequenceRC implements CharSequence, Comparable<CharSequence> {
    private final int lenLess1;
    private final CharSequence sequence;

    public SequenceRC( final CharSequence sequence ) {
        this.lenLess1 = sequence.length() - 1;
        this.sequence = sequence;
    }

    @Override
    public int length() {
        return sequence.length();
    }

    @Override
    public char charAt( final int index ) {
        return switch ( Character.toUpperCase(sequence.charAt(lenLess1 - index)) ) {
            case 'A' -> 'T';
            case 'B' -> 'V';
            case 'C' -> 'G';
            case 'D' -> 'H';
            case 'G' -> 'C';
            case 'H' -> 'D';
            case 'K' -> 'M';
            case 'M' -> 'K';
            case 'R' -> 'Y';
            case 'S' -> 'S';
            case 'T', 'U' -> 'A';
            case 'V' -> 'B';
            case 'W' -> 'W';
            case 'Y' -> 'R';
            default -> 'N';
        };
    }

    @Override
    public @NotNull CharSequence subSequence( final int start, final int end ) {
        return new StringBuilder(end - start).append(this, start, end).toString();
    }

    @Override
    public @NotNull String toString() {
        return new StringBuilder(this).toString();
    }

    @Override
    public int compareTo( final CharSequence charSequence ) {
        final int len1 = length();
        final int len2 = charSequence.length();
        final int cmpLen = Math.min(len1, len2);
        for ( int idx = 0; idx != cmpLen; ++idx ) {
            final char char1 = charAt(idx);
            final char char2 = Character.toUpperCase(charSequence.charAt(idx));
            if ( char1 > char2 ) return 1;
            if ( char1 < char2 ) return -1;
        }
        return Integer.compare(len1, len2);
    }
}
