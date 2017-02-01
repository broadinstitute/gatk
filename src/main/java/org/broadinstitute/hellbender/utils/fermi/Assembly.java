package org.broadinstitute.hellbender.utils.fermi;

import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;

import java.io.*;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;

/** an assembly is just a collection of contigs */
public final class Assembly {
    private final List<Contig> contigs;

    Assembly( final List<Contig> contigs ) { this.contigs = contigs; }

    public int getNContigs() { return contigs.size(); }
    public Contig getContig( final int idx ) { return contigs.get(idx); }
    public List<Contig> getContigs() { return Collections.unmodifiableList(contigs); }

    /** a sequence of bases, coverage data, and connections to other contigs */
    public final static class Contig {
        private final byte[] sequence;
        private final int nSupportingReads;
        private final byte[] perBaseCoverage;
        private List<Connection> connections;

        Contig( final byte[] sequence, final int nSupportingReads, final byte[] perBaseCoverage ) {
            this.sequence = sequence;
            this.nSupportingReads = nSupportingReads;
            this.perBaseCoverage = perBaseCoverage;
        }

        public byte[] getSequence() { return sequence; }
        public int getNSupportingReads() { return nSupportingReads; }
        public byte[] getPerBaseCoverage() { return perBaseCoverage; }
        public List<Connection> getConnections() { return Collections.unmodifiableList(connections); }

        void setConnections( final List<Connection> connections ) {
            this.connections = connections;
        }
    }

    /** a connection between contigs */
    public final static class Connection {
        private final Contig target;      // contig that overlaps the one that possesses this connection
        private final int overlapLen;     // bases in common -- negative overlap lengths are legal, and represent gaps
        private final boolean isRC;       // if target is a predecessor
        private final boolean isTargetRC; // if connection is to RC of target contig

        Connection( final Contig target, final int overlapLen, final boolean isRC, final boolean isTargetRC ) {
            this.target = target;
            this.overlapLen = overlapLen;
            this.isRC = isRC;
            this.isTargetRC = isTargetRC;
        }

        public Contig getTarget() { return target; }
        public int getOverlapLen() { return overlapLen; }
        public boolean isRC() { return isRC; }
        public boolean isTargetRC() { return isTargetRC; }
    }

    public void writeGFA( final String fileName, final PipelineOptions pipelineOptions ) {
        final HashMap<Contig, Integer> idMap = new HashMap<>((int)((contigs.size()*4L)/3) + 1);
        int id = 0;
        for (final Contig contig : contigs) {
            idMap.put(contig, ++id);
        }
        try ( final Writer writer =
                      new OutputStreamWriter(new BufferedOutputStream(BucketUtils.createFile(fileName, pipelineOptions))) ) {
            writer.write("H\tVN:Z:1\n");
            for ( Contig contig : contigs ) {
                final int contigId = idMap.get(contig);
                writer.write("S\ttig" + contigId + "\t" + new String(contig.getSequence()) +
                        "\tLN:i:" + contig.getSequence().length + "\tRC:i:" + contig.getNSupportingReads() + "\n");
                for ( Connection connection : contig.getConnections() ) {
                    final int targetId = idMap.get(connection.getTarget());
                    if ( contigId <= targetId ) {
                        final int overlapLen = connection.getOverlapLen();
                        writer.write("L\ttig" + contigId + "\t" + (connection.isRC() ? "-" : "+") +
                                "\ttig" + targetId + "\t" + (connection.isTargetRC() ? "-" : "+") + "\t" +
                                (overlapLen < 0 ? -overlapLen + "H" : overlapLen + "M") + "\n");
                    }
                }
            }
        } catch (final IOException ioe) {
            throw new GATKException("Can't write " + fileName, ioe);
        }
    }
}
