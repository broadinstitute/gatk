package org.broadinstitute.hellbender.utils.codecs;

import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.sv.DiscordantPairEvidence;
import org.broadinstitute.hellbender.utils.io.BlockCompressedIntervalStream.Reader;
import org.broadinstitute.hellbender.utils.io.BlockCompressedIntervalStream.Writer;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.IOException;
import java.util.List;

public class DiscordantPairEvidenceBCICodec extends AbstractBCICodec<DiscordantPairEvidence> {
    private boolean versionChecked = false;

    private static final String PE_BCI_FILE_EXTENSION = ".pe.bci";

    @Override
    public DiscordantPairEvidence decode( final Reader<DiscordantPairEvidence> reader )
            throws IOException {
        if ( !versionChecked ) {
            if ( !DiscordantPairEvidence.BCI_VERSION.equals(reader.getVersion()) ) {
                throw new UserException("baf.bci file has wrong version: expected " +
                        DiscordantPairEvidence.BCI_VERSION + " but found " + reader.getVersion());
            }
            versionChecked = true;
        }
        final DataInputStream dis = reader.getStream();
        final String sample = reader.getSampleNames().get(dis.readInt());
        final SAMSequenceDictionary dict = reader.getDictionary();
        final String startContig = dict.getSequence(dis.readInt()).getSequenceName();
        final int start = dis.readInt();
        final boolean startStrand = dis.readBoolean();
        final String endContig = dict.getSequence(dis.readInt()).getSequenceName();
        final int end = dis.readInt();
        final boolean endStrand = dis.readBoolean();
        return new DiscordantPairEvidence(sample, startContig, start, startStrand,
                                            endContig, end, endStrand);
    }

    @Override
    public Class<DiscordantPairEvidence> getFeatureType() { return DiscordantPairEvidence.class; }

    @Override
    public boolean canDecode( final String path ) {
        return path.toLowerCase().endsWith(PE_BCI_FILE_EXTENSION);
    }

    @Override
    public Writer<DiscordantPairEvidence> makeSink( final GATKPath path,
                                                    final SAMSequenceDictionary dict,
                                                    final List<String> sampleNames,
                                                    final int compressionLevel ) {
        final String className = DiscordantPairEvidence.class.getSimpleName();
        final FeaturesHeader header =
                new FeaturesHeader(className, DiscordantPairEvidence.BCI_VERSION, dict, sampleNames);
        return new Writer<>(path, header, this::encode, compressionLevel);
    }

    @Override
    public void encode( final DiscordantPairEvidence peEvidence,
                        final Writer<DiscordantPairEvidence> writer ) throws IOException {
        final DataOutputStream dos = writer.getStream();
        dos.writeInt(writer.getSampleIndex(peEvidence.getSample()));
        dos.writeInt(writer.getContigIndex(peEvidence.getContig()));
        dos.writeInt(peEvidence.getStart());
        dos.writeBoolean(peEvidence.getStartStrand());
        dos.writeInt(writer.getContigIndex(peEvidence.getEndContig()));
        dos.writeInt(peEvidence.getEndPosition());
        dos.writeBoolean(peEvidence.getEndStrand());
    }
}
