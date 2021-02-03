package org.broadinstitute.hellbender.utils.codecs;

import com.google.common.base.Splitter;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.IOUtil;
import htsjdk.tribble.AsciiFeatureCodec;
import htsjdk.tribble.index.tabix.TabixFormat;
import htsjdk.tribble.readers.LineIterator;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.tools.sv.LocusDepth;
import org.broadinstitute.hellbender.utils.io.FeatureOutputStream;

import java.util.List;

public class LocusDepthCodec extends AsciiFeatureCodec<LocusDepth>
        implements FeatureOutputCodec<LocusDepth, FeatureOutputStream<LocusDepth>> {
    public static final String FORMAT_SUFFIX = ".ld.txt";
    private static final Splitter splitter = Splitter.on("\t");

    public LocusDepthCodec() {
        super(LocusDepth.class);
    }

    @Override public TabixFormat getTabixFormat() {
        return new TabixFormat(TabixFormat.ZERO_BASED, 1, 2, 0, '#', 0);
    }

    @Override public LocusDepth decode( final String line ) {
        final List<String> tokens = splitter.splitToList(line);
        if ( tokens.size() != 7 ) {
            throw new IllegalArgumentException("Invalid number of columns: " + tokens.size());
        }
        return new LocusDepth(tokens.get(0),
                Integer.parseUnsignedInt(tokens.get(1)) + 1,
                (byte)tokens.get(2).charAt(0),
                Integer.parseUnsignedInt(tokens.get(3)),
                Integer.parseUnsignedInt(tokens.get(4)),
                Integer.parseUnsignedInt(tokens.get(5)),
                Integer.parseUnsignedInt(tokens.get(6)));
    }

    @Override public Object readActualHeader( LineIterator reader ) { return null; }

    @Override
    public boolean canDecode( final String path ) {
        String toDecode = path.toLowerCase();
        if ( IOUtil.hasBlockCompressedExtension(toDecode) ) {
            toDecode = toDecode.substring(0, toDecode.lastIndexOf('.'));
        }
        return toDecode.endsWith(FORMAT_SUFFIX);
    }

    @Override
    public FeatureOutputStream<LocusDepth> makeSink( final GATKPath path,
                                                     final SAMSequenceDictionary dict,
                                                     final List<String> sampleNames,
                                                     final int compressionLevel ) {
        return new FeatureOutputStream<>(path,
                getTabixFormat(),
                LocusDepthCodec::encode,
                dict,
                compressionLevel);
    }

    @Override
    public void encode( final LocusDepth ev, final FeatureOutputStream<LocusDepth> os ) {
        os.write(ev);
    }

    public static String encode( final LocusDepth locusDepth ) {
        return locusDepth.getContig() + "\t" + (locusDepth.getStart() - 1) + "\t" +
                locusDepth.getRefCall() + "\t" + locusDepth.getADepth() + "\t" +
                locusDepth.getCDepth() + "\t" + locusDepth.getGDepth() + "\t" +
                locusDepth.getTDepth();
    }
}
