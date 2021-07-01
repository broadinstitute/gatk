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
    private final int refIdx; // index into nucleotideValues
    private final int altIdx; // index into nucleotideValues
    private int totalDepth;
    private int altDepth;
    public final static String BCI_VERSION = "1.0";

    // our own private copy so that we don't make repeated array allocations
    private final static Nucleotide[] nucleotideValues = Nucleotide.values();

    public LocusDepth( final Locatable loc, final int refIdx, final int altIdx ) {
        this(loc.getContig(), loc.getStart(), refIdx, altIdx, 0, 0);
    }

    public LocusDepth( final String contig, final int position,
                       final int refIdx, final int altIdx,
                       final int totalDepth, final int altDepth ) {
        this.contig = contig;
        this.position = position;
        this.refIdx = refIdx;
        this.altIdx = altIdx;
        this.totalDepth = totalDepth;
        this.altDepth = altDepth;
    }

    public LocusDepth( final DataInputStream dis, final SAMSequenceDictionary dict ) throws IOException {
        contig = dict.getSequence(dis.readInt()).getSequenceName();
        position = dis.readInt();
        refIdx = dis.readByte();
        altIdx = dis.readByte();
        totalDepth = dis.readInt();
        altDepth = dis.readInt();
    }

    public void observe( final int idx ) {
        if ( idx == altIdx ) altDepth += 1;
        totalDepth += 1;
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
        return nucleotideValues[refIdx].encodeAsChar();
    }

    public char getAltCall() {
        return nucleotideValues[altIdx].encodeAsChar();
    }

    public int getAltDepth() {
        return altDepth;
    }

    public int getTotalDepth() {
        return totalDepth;
    }

    public void write( final DataOutputStream dos,
                       final SAMSequenceDictionary dict ) throws IOException {
        dos.writeInt(dict.getSequenceIndex(contig));
        dos.writeInt(position);
        dos.writeByte(refIdx);
        dos.writeByte(altIdx);
        dos.writeInt(totalDepth);
        dos.writeInt(altDepth);
    }

    public String toString() {
        return getContig() + "\t" + getStart() + "\t" + nucleotideValues[refIdx] + "\t" +
                nucleotideValues[altIdx] + "\t" + totalDepth + "\t" + altDepth;
    }
}
