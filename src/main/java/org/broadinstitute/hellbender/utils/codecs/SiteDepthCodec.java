package org.broadinstitute.hellbender.utils.codecs;

import com.google.common.base.Splitter;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.IOUtil;
import htsjdk.tribble.AsciiFeatureCodec;
import htsjdk.tribble.index.tabix.TabixFormat;
import htsjdk.tribble.readers.LineIterator;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.sv.SiteDepth;
import org.broadinstitute.hellbender.tools.sv.SiteDepthSortMerger;
import org.broadinstitute.hellbender.utils.io.FeatureOutputStream;

import java.util.List;

/** Codec to handle SiteDepths in tab-delimited text files */
public class SiteDepthCodec extends AsciiFeatureCodec<SiteDepth>
        implements FeatureOutputCodec<SiteDepth, FeatureOutputStream<SiteDepth>> {
    public static final String FORMAT_SUFFIX = ".sd.txt";
    private static final Splitter splitter = Splitter.on("\t");

    public SiteDepthCodec() {
        super(SiteDepth.class);
    }

    @Override public TabixFormat getTabixFormat() {
        return new TabixFormat(TabixFormat.ZERO_BASED, 1, 2, 0, '#', 0);
    }

    @Override public SiteDepth decode( final String line ) {
        final List<String> tokens = splitter.splitToList(line);
        if ( tokens.size() != 7 ) {
            throw new IllegalArgumentException("Invalid number of columns: " + tokens.size());
        }
        return new SiteDepth(tokens.get(0),
                Integer.parseUnsignedInt(tokens.get(1)) + 1,
                tokens.get(2),
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
    public FeatureOutputStream<SiteDepth> makeSink( final GATKPath path,
                                                    final SAMSequenceDictionary dict,
                                                    final List<String> sampleNames,
                                                    final int compressionLevel ) {
        if ( sampleNames.size() != 1 ) {
            throw new UserException("SiteDepth records do not encode their sample, and must all " +
                    "refer to a single sample, but the list of sample names is of " +
                    "size=" + sampleNames.size());
        }
        return new FeatureOutputStream<>(path,
                getTabixFormat(),
                SiteDepthCodec::encode,
                dict,
                compressionLevel);
    }

    @Override
    public void encode( final SiteDepth ev, final FeatureOutputStream<SiteDepth> os ) {
        os.write(ev);
    }

    @Override
    public FeatureSink<SiteDepth> makeSortMerger( final GATKPath path,
                                                  final SAMSequenceDictionary dict,
                                                  final List<String> sampleNames,
                                                  final int compressionLevel ) {
        return new SiteDepthSortMerger(dict, makeSink(path, dict, sampleNames, compressionLevel));
    }

    public static String encode( final SiteDepth siteDepth ) {
        return siteDepth.getContig() + "\t" + (siteDepth.getStart() - 1) + "\t" +
                siteDepth.getSample() + "\t" +
                siteDepth.getDepth(0) + "\t" + siteDepth.getDepth(1) + "\t" +
                siteDepth.getDepth(2) + "\t" + siteDepth.getDepth(3);
    }
}
