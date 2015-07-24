package org.broadinstitute.hellbender.utils.read.markduplicates;

import htsjdk.samtools.util.SortingCollection;
import org.broadinstitute.hellbender.exceptions.GATKException;


import java.io.*;

/** Coded for ReadEnds that just outputs the primitive fields and reads them back. */
public final class ReadEndsForMarkDuplicatesCodec implements SortingCollection.Codec<ReadEndsForMarkDuplicates> {
    private DataInputStream in;
    private DataOutputStream out;

    /**
     * For an explanation of why codecs must implement clone(),
     * see the HTSJDK documentation for {@link SortingCollection.Codec}.
     */
    @Override
    public ReadEndsForMarkDuplicatesCodec clone() {
        return new ReadEndsForMarkDuplicatesCodec();
    }

    public void setOutputStream(final OutputStream os) { this.out = new DataOutputStream(os); }

    public void setInputStream(final InputStream is) { this.in = new DataInputStream(is); }

    public DataInputStream getInputStream() {
        return in;
    }

    public DataOutputStream getOutputStream() {
        return out;
    }

    public void encode(final ReadEndsForMarkDuplicates read) {
        try {
            this.out.writeShort(read.score);
            this.out.writeShort(read.libraryId);
            this.out.writeByte(read.orientation);
            this.out.writeInt(read.read1ReferenceIndex);
            this.out.writeInt(read.read1Coordinate);
            this.out.writeLong(read.read1IndexInFile);
            this.out.writeInt(read.read2ReferenceIndex);

            if (read.orientation > ReadEnds.R) {
                this.out.writeInt(read.read2Coordinate);
                this.out.writeLong(read.read2IndexInFile);
            }

            this.out.writeShort(read.readGroup);
            this.out.writeShort(read.tile);
            this.out.writeShort(read.x);
            this.out.writeShort(read.y);
            this.out.writeByte(read.orientationForOpticalDuplicates);
        } catch (final IOException ioe) {
            throw new GATKException("Exception writing ReadEnds to file.", ioe);
        }
    }

    public ReadEndsForMarkDuplicates decode() {
        final ReadEndsForMarkDuplicates read = new ReadEndsForMarkDuplicates();
        try {
            // If the first read results in an EOF we've exhausted the stream
            try {
                read.score = this.in.readShort();
            } catch (final EOFException eof) {
                return null;
            }

            read.libraryId = this.in.readShort();
            read.orientation = this.in.readByte();
            read.read1ReferenceIndex = this.in.readInt();
            read.read1Coordinate = this.in.readInt();
            read.read1IndexInFile = this.in.readLong();
            read.read2ReferenceIndex = this.in.readInt();

            if (read.orientation > ReadEnds.R) {
                read.read2Coordinate = this.in.readInt();
                read.read2IndexInFile = this.in.readLong();
            }

            read.readGroup = this.in.readShort();
            read.tile = this.in.readShort();
            read.x = this.in.readShort();
            read.y = this.in.readShort();

            read.orientationForOpticalDuplicates = this.in.readByte();

            return read;
        } catch (final IOException ioe) {
            throw new GATKException("Exception writing ReadEnds to file.", ioe);
        }
    }
}
