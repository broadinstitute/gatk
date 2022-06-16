package org.broadinstitute.hellbender.utils.codecs;

import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.sv.BafEvidence;
import org.broadinstitute.hellbender.tools.sv.BafEvidenceSortMerger;
import org.broadinstitute.hellbender.tools.sv.SVFeaturesHeader;
import org.broadinstitute.hellbender.utils.io.BlockCompressedIntervalStream.Reader;
import org.broadinstitute.hellbender.utils.io.BlockCompressedIntervalStream.Writer;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.IOException;
import java.util.List;

/** Codec to handle BafEvidence in BlockCompressedInterval files */
public class BafEvidenceBCICodec extends AbstractBCICodec<BafEvidence> {
    private boolean versionChecked = false;
    private static final String BAF_BCI_FILE_EXTENSION = ".baf.bci";

    @Override
    public BafEvidence decode( final Reader<BafEvidence> reader ) throws IOException {
        if ( !versionChecked ) {
            if ( !BafEvidence.BCI_VERSION.equals(reader.getVersion()) ) {
                throw new UserException("baf.bci file has wrong version: expected " +
                        BafEvidence.BCI_VERSION + " but found " + reader.getVersion());
            }
            versionChecked = true;
        }
        final DataInputStream dis = reader.getStream();
        final String sample = reader.getSampleNames().get(dis.readInt());
        final String contig = reader.getDictionary().getSequence(dis.readInt()).getSequenceName();
        final int position = dis.readInt();
        final double value = dis.readDouble();
        return new BafEvidence(sample, contig, position, value);
    }

    @Override
    public Class<BafEvidence> getFeatureType() { return BafEvidence.class; }

    @Override
    public boolean canDecode( final String path ) {
        return path.toLowerCase().endsWith(BAF_BCI_FILE_EXTENSION);
    }

    @Override
    public Writer<BafEvidence> makeSink( final GATKPath path,
                                         final SAMSequenceDictionary dict,
                                         final List<String> sampleNames,
                                         final int compressionLevel ) {
        final String className = BafEvidence.class.getSimpleName();
        return new Writer<>(path,
                            new SVFeaturesHeader(className, BafEvidence.BCI_VERSION, dict, sampleNames),
                            this::encode,
                            compressionLevel);
    }

    @Override
    public void encode( final BafEvidence bafEvidence, final Writer<BafEvidence> writer ) throws IOException {
        final DataOutputStream dos = writer.getStream();
        dos.writeInt(writer.getSampleIndex(bafEvidence.getSample()));
        dos.writeInt(writer.getContigIndex(bafEvidence.getContig()));
        dos.writeInt(bafEvidence.getStart());
        dos.writeDouble(bafEvidence.getValue());
    }

    @Override
    public FeatureSink<BafEvidence> makeSortMerger( final GATKPath path,
                                                    final SAMSequenceDictionary dict,
                                                    final List<String> sampleNames,
                                                    final int compressionLevel ) {
        return new BafEvidenceSortMerger(dict, makeSink(path, dict, sampleNames, compressionLevel));
    }
}
