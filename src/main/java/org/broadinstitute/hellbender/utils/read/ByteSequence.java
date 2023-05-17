package org.broadinstitute.hellbender.utils.read;

import org.jetbrains.annotations.NotNull;

/**
 * Immutable sequence of bytes (assuming you don't change the array externally).
 * Subsequences reference the parent array.
 */
public final class ByteSequence implements CharSequence {
    public static final ByteSequence EMPTY = new ByteSequence(new byte[0]);

    private final byte[] bytes;
    private final int start;
    private final int length;

    public ByteSequence( final byte[] bytes ) {
        this.bytes = bytes;
        this.start = 0;
        this.length = bytes.length;
    }

    public ByteSequence( final ByteSequence bSeq, final int start, final int length ) {
        if ( start < 0 || length < 0 || start + length > bSeq.length ) {
            throw new IndexOutOfBoundsException();
        }
        this.bytes = bSeq.bytes;
        this.start = bSeq.start + start;
        this.length = length;
    }

    public ByteSequence( final byte val ) {
        this.bytes = new byte[1];
        bytes[0] = val;
        this.start = 0;
        this.length = 1;
    }

    public int getStart() { return start; }

    @Override
    public int length() { return length; }

    public byte byteAt( final int idx ) {
        if ( idx < 0 || idx >= length ) {
            throw new IndexOutOfBoundsException();
        }
        return bytes[start + idx];
    }

    @Override
    public char charAt( final int idx ) {
        return (char)byteAt(idx);
    }

    public ByteSequence subSequence( final int start ) {
        return new ByteSequence(this, start, length() - start);
    }

    /**
     * Note: args are start+length, not start+end
    public @NotNull ByteSequence subSequence( final int start, final int len ) {
        return new ByteSequence(this, start, len);
    }
     */

    public ByteSequence append( final ByteSequence bytesToAppend ) {
        final byte[] newBytes = new byte[length + bytesToAppend.length];
        System.arraycopy(bytes, start, newBytes, 0, length);
        System.arraycopy(bytesToAppend.bytes, bytesToAppend.start, newBytes, length, bytesToAppend.length);
        return new ByteSequence(newBytes);
    }

    @Override
    public boolean equals( final Object obj ) {
        if ( !(obj instanceof final ByteSequence that) ) {
            return false;
        }
        if ( this.length != that.length ) {
            return false;
        }
        for ( int idx = 0; idx != length; ++idx ) {
            if ( this.byteAt(idx) != that.byteAt(idx) ) {
                return false;
            }
        }
        return true;
    }

    @Override
    public int hashCode() {
        int hash = 83;
        for ( int idx = 0; idx != length; ++idx ) {
            hash = (47 * hash) + bytes[idx];
        }
        return 41 * hash;
    }

    @Override
    public @NotNull String toString() {
        return new String(bytes, start, length);
    }
}
