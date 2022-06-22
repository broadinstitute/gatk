package org.broadinstitute.hellbender.utils.codecs;

import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.sv.SiteDepth;
import org.broadinstitute.hellbender.tools.sv.SiteDepthSortMerger;
import org.broadinstitute.hellbender.tools.sv.SVFeaturesHeader;
import org.broadinstitute.hellbender.utils.io.BlockCompressedIntervalStream.Reader;
import org.broadinstitute.hellbender.utils.io.BlockCompressedIntervalStream.Writer;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.IOException;
import java.util.List;

/** Codec to handle SiteDepths in BlockCompressedInterval files */
public class SiteDepthBCICodec extends AbstractBCICodec<SiteDepth> {
    private boolean versionChecked = false;
    private static final String SD_BCI_FILE_EXTENSION = ".sd.bci";

    @Override
    public SiteDepth decode( final Reader<SiteDepth> reader ) throws IOException {
        if ( !versionChecked ) {
            if ( !SiteDepth.BCI_VERSION.equals(reader.getVersion()) ) {
                throw new UserException("bci file has wrong version: expected " +
                        SiteDepth.BCI_VERSION + " but found " + reader.getVersion());
            }
            versionChecked = true;
        }
        final DataInputStream dis = reader.getStream();
        return new SiteDepth(reader.getDictionary().getSequence(dis.readInt()).getSequenceName(),
                                dis.readInt(), reader.getSampleNames().get(dis.readInt()),
                                dis.readInt(), dis.readInt(), dis.readInt(), dis.readInt());
    }

    @Override
    public Class<SiteDepth> getFeatureType() { return SiteDepth.class; }

    @Override
    public boolean canDecode( final String path ) {
        return path.toLowerCase().endsWith(SD_BCI_FILE_EXTENSION);
    }

    @Override
    public Writer<SiteDepth> makeSink( final GATKPath path,
                                       final SAMSequenceDictionary dict,
                                       final List<String> sampleNames,
                                       final int compressionLevel ) {
        final String className = SiteDepth.class.getSimpleName();
        return new Writer<>(path,
                            new SVFeaturesHeader(className, SiteDepth.BCI_VERSION, dict, sampleNames),
                            this::encode,
                            compressionLevel);
    }

    @Override
    public void encode( final SiteDepth siteDepth, final Writer<SiteDepth> writer )
            throws IOException {
        final DataOutputStream dos = writer.getStream();
        dos.writeInt(writer.getContigIndex(siteDepth.getContig()));
        dos.writeInt(siteDepth.getStart());
        dos.writeInt(writer.getSampleIndex(siteDepth.getSample()));
        dos.writeInt(siteDepth.getDepth(0));
        dos.writeInt(siteDepth.getDepth(1));
        dos.writeInt(siteDepth.getDepth(2));
        dos.writeInt(siteDepth.getDepth(3));
    }

    @Override
    public FeatureSink<SiteDepth> makeSortMerger( final GATKPath path,
                                                  final SAMSequenceDictionary dict,
                                                  final List<String> sampleNames,
                                                  final int compressionLevel ) {
        return new SiteDepthSortMerger(dict, makeSink(path, dict, sampleNames, compressionLevel));
    }
}
