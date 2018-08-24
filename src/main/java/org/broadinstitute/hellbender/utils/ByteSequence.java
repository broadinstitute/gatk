package org.broadinstitute.hellbender.utils;

import java.util.Iterator;
import java.util.NoSuchElementException;

public interface ByteSequence extends Iterable<Byte> {
    byte byteAt( int index );
    int length();
    ByteSequence subSequence( int start, int end );

    @Override default ByteIterator iterator() { return iterator(0); }
    default ByteIterator iterator( final int index ) { return new ByteSequenceIterator(this, index); }
    default ByteIterator reverseIterator() { return reverseIterator(0); }
    default ByteIterator reverseIterator( final int index ) { return new ByteSequenceReverseIterator(this, index); }
    default ByteIterator rcIterator( final Complementer complementer ) { return rcIterator(0, complementer); }
    default ByteIterator rcIterator( final int index, final Complementer complementer ) {
        return new ByteSequenceReverseComplementIterator(this, index, complementer);
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
            this(byteSequence, 0);
        }

        public ByteSequenceIterator( final ByteSequence byteSequence, final int index ) {
            Utils.nonNull(byteSequence);
            Utils.validateArg(index >= 0 && index < byteSequence.length(), "index out of bounds");
            this.byteSequence = byteSequence;
            this.index = index;
        }

        @Override public boolean hasNext() { return index < byteSequence.length(); }
        @Override public byte nextByte() { return byteSequence.byteAt(index++); }
    }

    final class ByteSequenceReverseIterator implements ByteIterator {
        private final ByteSequence byteSequence;
        private int index;

        public ByteSequenceReverseIterator( final ByteSequence byteSequence ) {
            this(byteSequence, 0);
        }

        public ByteSequenceReverseIterator( final ByteSequence byteSequence, final int index ) {
            Utils.nonNull(byteSequence);
            Utils.validateArg(index >= 0 && index < byteSequence.length(), "index out of bounds");
            this.byteSequence = byteSequence;
            this.index = byteSequence.length() - index;
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
            this(byteSequence, 0, complementer);
        }

        public ByteSequenceReverseComplementIterator( final ByteSequence byteSequence,
                                                      final int index,
                                                      final Complementer complementer ) {
            Utils.nonNull(byteSequence);
            Utils.validateArg(index >= 0 && index < byteSequence.length(), "index out of bounds");
            this.byteSequence = byteSequence;
            this.index = byteSequence.length() - index;
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
