package org.broadinstitute.hellbender.utils.codecs;

import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.sv.LocusDepth;
import org.broadinstitute.hellbender.tools.sv.LocusDepthSortMerger;
import org.broadinstitute.hellbender.tools.sv.SVFeaturesHeader;
import org.broadinstitute.hellbender.utils.io.BlockCompressedIntervalStream.Reader;
import org.broadinstitute.hellbender.utils.io.BlockCompressedIntervalStream.Writer;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.IOException;
import java.util.List;

/** Codec to handle LocusDepths in BlockCompressedInterval files */
public class LocusDepthBCICodec extends AbstractBCICodec<LocusDepth> {
    private boolean versionChecked = false;
    private static final String LD_BCI_FILE_EXTENSION = ".ld.bci";

    @Override
    public LocusDepth decode( final Reader<LocusDepth> reader ) throws IOException {
        if ( !versionChecked ) {
            if ( !LocusDepth.BCI_VERSION.equals(reader.getVersion()) ) {
                throw new UserException("bci file has wrong version: expected " +
                        LocusDepth.BCI_VERSION + " but found " + reader.getVersion());
            }
            versionChecked = true;
        }
        final DataInputStream dis = reader.getStream();
        return new LocusDepth(reader.getDictionary().getSequence(dis.readInt()).getSequenceName(),
                                dis.readInt(),
                                reader.getSampleNames().get(dis.readInt()),
                                dis.readByte(), dis.readByte(),
                                dis.readInt(), dis.readInt(), dis.readInt(), dis.readInt());
    }

    @Override
    public Class<LocusDepth> getFeatureType() { return LocusDepth.class; }

    @Override
    public boolean canDecode( final String path ) {
        return path.toLowerCase().endsWith(LD_BCI_FILE_EXTENSION);
    }

    @Override
    public Writer<LocusDepth> makeSink( final GATKPath path,
                                        final SAMSequenceDictionary dict,
                                        final List<String> sampleNames,
                                        final int compressionLevel ) {
        if ( sampleNames.size() != 1 ) {
            throw new UserException("LocusDepth records do not encode their sample, and must all " +
                                    "refer to a single sample, but the list of sample names is of " +
                                    "size=" + sampleNames.size());
        }
        final String className = LocusDepth.class.getSimpleName();
        return new Writer<>(path,
                            new SVFeaturesHeader(className, LocusDepth.BCI_VERSION, dict, sampleNames),
                            this::encode,
                            compressionLevel);
    }

    @Override
    public void encode( final LocusDepth locusDepth, final Writer<LocusDepth> writer )
            throws IOException {
        final DataOutputStream dos = writer.getStream();
        dos.writeInt(writer.getContigIndex(locusDepth.getContig()));
        dos.writeInt(locusDepth.getStart());
        dos.writeInt(writer.getSampleIndex(locusDepth.getSample()));
        dos.writeByte(locusDepth.getRefIndex());
        dos.writeByte(locusDepth.getAltIndex());
        dos.writeInt(locusDepth.getDepth(0));
        dos.writeInt(locusDepth.getDepth(1));
        dos.writeInt(locusDepth.getDepth(2));
        dos.writeInt(locusDepth.getDepth(3));
    }

    @Override
    public FeatureSink<LocusDepth> makeSortMerger( final GATKPath path,
                                                   final SAMSequenceDictionary dict,
                                                   final List<String> sampleNames,
                                                   final int compressionLevel ) {
        return new LocusDepthSortMerger(dict, makeSink(path, dict, sampleNames, compressionLevel));
    }
}
