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

}
