package org.broadinstitute.hellbender.utils.read;

import htsjdk.samtools.SAMUtils;

import java.io.BufferedWriter;
import java.io.IOException;

public final class UnalignedRead {
    private final String name;
    private final ByteSequence calls;
    private final ByteSequence quals;

    public UnalignedRead( final String name, final ByteSequence calls, final ByteSequence quals ) {
        this.name = name;
        this.calls = calls;
        this.quals = quals;
    }

    public UnalignedRead( final GATKRead read ) {
        name = read.getName();
        calls = new ByteSequence(read.getBasesNoCopy());
        quals = new ByteSequence(read.getBaseQualitiesNoCopy());
    }

    public String getName() { return name; }
    public ByteSequence getCalls() { return calls; }
    public ByteSequence getQuals() { return quals; }

    public void writeFASTQ( final BufferedWriter writer ) throws IOException {
        writer.write("@");
        writer.write(name);
        writer.newLine();
        int lineLength = 0;
        final int length = calls.length();
        for ( int idx = 0; idx != length; ++idx ) {
            writer.write(calls.byteAt(idx));
            if ( ++lineLength == 80 ) {
                writer.newLine();
                lineLength = 0;
            }
        }
        if ( lineLength > 0 ) {
            writer.newLine();
            lineLength = 0;
        }
        writer.write("+");
        writer.newLine();
        for ( int idx = 0; idx != length; ++idx ) {
            writer.write(SAMUtils.phredToFastq(quals.byteAt(idx)));
            if ( ++lineLength == 80 ) {
                writer.newLine();
                lineLength = 0;
            }
        }
        if ( lineLength > 0 ) {
            writer.newLine();
        }
    }

    public void writeFASTA( final BufferedWriter writer ) throws IOException {
        writer.write(">");
        writer.write(name);
        writer.newLine();
        int lineLength = 0;
        final int length = calls.length();
        for ( int idx = 0; idx != length; ++idx ) {
            writer.write(calls.byteAt(idx));
            if ( ++lineLength == 80 ) {
                writer.newLine();
                lineLength = 0;
            }
        }
        if ( lineLength > 0 ) {
            writer.newLine();
        }
    }

    /**
     * Immutable sequence of bytes (assuming you don't change the array directly).
     * Subsequences reference the parent array.
     */
    public static final class ByteSequence {
        public static final ByteSequence EMPTY = new ByteSequence(new byte[0]);
        final byte[] bytes;
        final int start;
        final int length;

        public ByteSequence( final byte[] bytes ) {
            this.bytes = bytes;
            this.start = 0;
            this.length = bytes.length;
        }

        private ByteSequence( final ByteSequence bSeq, final int start, final int length ) {
            if ( start < 0 || length < 0 || start + length > bSeq.length ) {
                throw new IndexOutOfBoundsException();
            }
            this.bytes = bSeq.bytes;
            this.start = bSeq.start + start;
            this.length = length;
        }

        public int getStart() { return start; }
        public int length() { return length; }

        public byte byteAt( final int idx ) {
            if ( idx < 0 || idx >= length ) {
                throw new IndexOutOfBoundsException();
            }
            return bytes[start + idx];
        }

        public ByteSequence subSequence( final int start ) {
            return new ByteSequence(this, start, length() - start);
        }

        /** Note: args are start+length, not start+end */
        public ByteSequence subSequence( final int start, final int len ) {
            return new ByteSequence( this, start, len);
        }

        @Override
        public boolean equals( final Object obj ) {
            if ( !(obj instanceof ByteSequence) ) {
                return false;
            }
            final ByteSequence that = (ByteSequence)obj;
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
        public String toString() { return new String(bytes, start, length); }
    }
}
