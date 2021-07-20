package org.broadinstitute.hellbender.tools.sv;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Locatable;
import htsjdk.tribble.Feature;
import org.broadinstitute.hellbender.utils.Nucleotide;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.IOException;

@VisibleForTesting
public final class LocusDepth implements Feature {
    private final String contig;
    private final int position;
    private final byte refCall; // index into nucleotideValues
    private int[] depths;
    public final static String BCI_VERSION = "1.0";

    // our own private copy so that we don't make repeated array allocations
    private final static Nucleotide[] nucleotideValues = Nucleotide.values();

    public LocusDepth( final Locatable loc, final byte refCall ) {
        this.contig = loc.getContig();
        this.position = loc.getStart();
        this.refCall = refCall;
        this.depths = new int[4];
    }

    public LocusDepth( final String contig, final int position, final byte refCall,
                       final int aDepth, final int cDepth, final int gDepth, final int tDepth ) {
        this.contig = contig;
        this.position = position;
        this.refCall = refCall;
        this.depths = new int[4];
        depths[0] = aDepth;
        depths[1] = cDepth;
        depths[2] = gDepth;
        depths[3] = tDepth;
    }

    public LocusDepth( final DataInputStream dis, final SAMSequenceDictionary dict ) throws IOException {
        contig = dict.getSequence(dis.readInt()).getSequenceName();
        position = dis.readInt();
        refCall = dis.readByte();
        depths = new int[4];
        depths[0] = dis.readInt();
        depths[1] = dis.readInt();
        depths[2] = dis.readInt();
        depths[3] = dis.readInt();
    }

    public void observe( final int idx ) {
        depths[idx] += 1;
    }

    @Override
    public String getContig() {
        return contig;
    }

    @Override
    public int getEnd() {
        return position;
    }

    @Override
    public int getStart() {
        return position;
    }

    public char getRefCall() {
        return (char)refCall;
    }

    public int getADepth() { return depths[0]; }
    public int getCDepth() { return depths[1]; }
    public int getGDepth() { return depths[2]; }
    public int getTDepth() { return depths[3]; }

    public void write( final DataOutputStream dos,
                       final SAMSequenceDictionary dict ) throws IOException {
        dos.writeInt(dict.getSequenceIndex(contig));
        dos.writeInt(position);
        dos.writeByte(refCall);
        dos.writeInt(depths[0]);
        dos.writeInt(depths[1]);
        dos.writeInt(depths[2]);
        dos.writeInt(depths[3]);
    }

    public String toString() {
        return contig + "\t" + position + "\t" + (char)refCall + "\t" +
                depths[0] + "\t" + depths[1] + "\t" + depths[2] + "\t" + depths[3];
    }
}
