package org.broadinstitute.hellbender.utils.codecs;

import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.sv.DepthEvidence;
import org.broadinstitute.hellbender.utils.io.BlockCompressedIntervalStream.Reader;
import org.broadinstitute.hellbender.utils.io.BlockCompressedIntervalStream.Writer;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.IOException;
import java.util.List;

public class DepthEvidenceBCICodec extends AbstractBCICodec<DepthEvidence> {
    private boolean versionChecked = false;
    private static final String RD_BCI_FILE_EXTENSION = ".rd.bci";

    @Override
    public DepthEvidence decode( final Reader<DepthEvidence> reader ) throws IOException {
        if ( !versionChecked ) {
            if ( !DepthEvidence.BCI_VERSION.equals(reader.getVersion()) ) {
                throw new UserException("baf.bci file has wrong version: expected " +
                        DepthEvidence.BCI_VERSION + " but found " + reader.getVersion());
            }
            versionChecked = true;
        }
        final DataInputStream dis = reader.getStream();
        final String contig = reader.getDictionary().getSequence(dis.readInt()).getSequenceName();
        final int start = dis.readInt();
        final int end = dis.readInt();
        final int nCounts = dis.readInt();
        final int[] counts = new int[nCounts];
        for ( int idx = 0; idx != nCounts; ++idx ) {
            counts[idx] = dis.readInt();
        }
        return new DepthEvidence(contig, start, end, counts);
    }

    @Override
    public Class<DepthEvidence> getFeatureType() { return DepthEvidence.class; }

    @Override
    public boolean canDecode( final String path ) {
        return path.toLowerCase().endsWith(RD_BCI_FILE_EXTENSION);
    }

    @Override
    public Writer<DepthEvidence> makeSink( final GATKPath path,
                                           final SAMSequenceDictionary dict,
                                           final List<String> sampleNames,
                                           final int compressionLevel ) {
        final String className = DepthEvidence.class.getSimpleName();
        return new Writer<>(path,
                            new FeaturesHeader(className, DepthEvidence.BCI_VERSION, dict, sampleNames),
                            this::encode,
                            compressionLevel);
    }

    @Override
    public void encode( final DepthEvidence depthEvidence,
                        final Writer<DepthEvidence> writer ) throws IOException {
        final DataOutputStream dos = writer.getStream();
        dos.writeInt(writer.getContigIndex(depthEvidence.getContig()));
        dos.writeInt(depthEvidence.getStart());
        dos.writeInt(depthEvidence.getEnd());
        final int[] counts = depthEvidence.getCounts();
        dos.writeInt(counts.length);
        for ( final int count : counts ) {
            dos.writeInt(count);
        }
    }
}
