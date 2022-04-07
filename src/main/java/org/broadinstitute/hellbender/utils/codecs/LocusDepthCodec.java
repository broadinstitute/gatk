package org.broadinstitute.hellbender.utils.codecs;

import com.google.common.base.Splitter;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.IOUtil;
import htsjdk.tribble.AsciiFeatureCodec;
import htsjdk.tribble.index.tabix.TabixFormat;
import htsjdk.tribble.readers.LineIterator;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.sv.LocusDepth;
import org.broadinstitute.hellbender.tools.sv.LocusDepthSortMerger;
import org.broadinstitute.hellbender.utils.io.FeatureOutputStream;

import java.util.List;

/** Codec to handle LocusDepths in tab-delimited text files */
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
        if ( tokens.size() != 8 ) {
            throw new IllegalArgumentException("Invalid number of columns: " + tokens.size());
        }
        return new LocusDepth(tokens.get(0),
                Integer.parseUnsignedInt(tokens.get(1)) + 1,
                tokens.get(2),
                (byte)tokens.get(3).charAt(0),
                Integer.parseUnsignedInt(tokens.get(4)),
                Integer.parseUnsignedInt(tokens.get(5)),
                Integer.parseUnsignedInt(tokens.get(6)),
                Integer.parseUnsignedInt(tokens.get(7)));
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
        if ( sampleNames.size() != 1 ) {
            throw new UserException("LocusDepth records do not encode their sample, and must all " +
                    "refer to a single sample, but the list of sample names is of " +
                    "size=" + sampleNames.size());
        }
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

    @Override
    public FeatureSink<LocusDepth> makeSortMerger( final GATKPath path,
                                                   final SAMSequenceDictionary dict,
                                                   final List<String> sampleNames,
                                                   final int compressionLevel ) {
        return new LocusDepthSortMerger(dict, makeSink(path, dict, sampleNames, compressionLevel));
    }

    public static String encode( final LocusDepth locusDepth ) {
        return locusDepth.getContig() + "\t" + (locusDepth.getStart() - 1) + "\t" +
                locusDepth.getSample() + "\t" + locusDepth.getRefCall() + "\t" +
                locusDepth.getADepth() + "\t" + locusDepth.getCDepth() + "\t" +
                locusDepth.getGDepth() + "\t" + locusDepth.getTDepth();
    }
}
