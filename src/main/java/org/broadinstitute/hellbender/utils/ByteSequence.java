package org.broadinstitute.hellbender.utils;

import java.util.Iterator;
import java.util.NoSuchElementException;

public interface ByteSequence extends Iterable<Byte> {
    byte byteAt( int index );
    int length();
    ByteSequence subSequence( int start, int end );

    @Override default ByteIterator iterator() { return new ByteSequenceIterator(this); }
    default ByteIterator reverseIterator() { return new ByteSequenceReverseIterator(this); }
    default ByteIterator rcIterator( final Complementer complementer ) {
        return new ByteSequenceReverseComplementIterator(this, complementer);
    }

    interface ByteIterator extends Iterator<Byte> {
        @Override default Byte next() {
            if ( !hasNext() ) {
                throw new NoSuchElementException("ByteIterator exhausted");
            }
            return nextByte();
        }
        byte nextByte();
    }

    interface Complementer {
        byte complement( byte val );
    }

    final class ByteSequenceIterator implements ByteIterator {
        private final ByteSequence byteSequence;
        private int index;

        public ByteSequenceIterator( final ByteSequence byteSequence ) {
            Utils.nonNull(byteSequence);
            this.byteSequence = byteSequence;
            this.index = 0;
        }

        @Override public boolean hasNext() { return index < byteSequence.length(); }
        @Override public byte nextByte() { return byteSequence.byteAt(index++); }
    }

    final class ByteSequenceReverseIterator implements ByteIterator {
        private final ByteSequence byteSequence;
        private int index;

        public ByteSequenceReverseIterator( final ByteSequence byteSequence ) {
            Utils.nonNull(byteSequence);
            this.byteSequence = byteSequence;
            this.index = byteSequence.length();
        }

        @Override public boolean hasNext() { return index > 0; }
        @Override public byte nextByte() { return byteSequence.byteAt(--index); }
    }

    final class ByteSequenceReverseComplementIterator implements ByteIterator {
        private final ByteSequence byteSequence;
        private int index;
        private final Complementer complementer;

        public ByteSequenceReverseComplementIterator( final ByteSequence byteSequence,
                                                      final Complementer complementer ) {
            Utils.nonNull(byteSequence);
            this.byteSequence = byteSequence;
            this.index = byteSequence.length();
            this.complementer = complementer;
        }

        @Override public boolean hasNext() { return index > 0; }
        @Override public byte nextByte() { return complementer.complement(byteSequence.byteAt(--index)); }
    }

    final class ByteArraySequence implements ByteSequence {
        private final byte[] bytes;

        public ByteArraySequence( final byte[] bytes ) {
            Utils.nonNull(bytes);
            this.bytes = bytes;
        }

        @Override public byte byteAt( final int index ) { return bytes[index]; }

        @Override public int length() { return bytes.length; }

        @Override public ByteSequence subSequence( final int start, final int end ) {
            return new ByteArraySubSequence(bytes, start, end);
        }
    }

    final class ByteArraySubSequence implements ByteSequence {
        private final byte[] bytes;
        private final int start;
        private final int length;

        public ByteArraySubSequence( final byte[] bytes, final int start, final int end ) {
            Utils.nonNull(bytes);
            Utils.validateArg(start >= 0, "start must be non-negative");
            Utils.validateArg(end <= bytes.length, "end must be no greater than array length");
            Utils.validateArg(end >= start, "end must be at least as large as start");
            this.bytes = bytes;
            this.start = start;
            this.length = end - start;
        }

        @Override public byte byteAt( final int index ) {
            Utils.validateArg(index >= 0 && index < length, "index out of bounds");
            return bytes[start + index];
        }

        @Override public int length() { return length; }

        @Override public ByteSequence subSequence( final int start, final int end ) {
            return new ByteArraySubSequence(bytes, this.start + start, this.start + end);
        }
    }
}
