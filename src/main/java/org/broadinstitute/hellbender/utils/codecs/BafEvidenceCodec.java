package org.broadinstitute.hellbender.utils.codecs;

import com.google.common.base.Splitter;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.IOUtil;
import htsjdk.tribble.AsciiFeatureCodec;
import htsjdk.tribble.index.tabix.TabixFormat;
import htsjdk.tribble.readers.LineIterator;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.tools.sv.BafEvidence;
import org.broadinstitute.hellbender.utils.io.FeatureOutputStream;

import java.util.Arrays;
import java.util.List;

public class BafEvidenceCodec extends AsciiFeatureCodec<BafEvidence>
        implements FeatureOutputCodec<BafEvidence, FeatureOutputStream<BafEvidence>> {

    public static final String FORMAT_SUFFIX = ".baf.txt";
    public static final String COL_DELIMITER = "\t";
    private static final Splitter splitter = Splitter.on(COL_DELIMITER);

    public BafEvidenceCodec() {
        super(BafEvidence.class);
    }

    @Override
    public TabixFormat getTabixFormat() {
        return new TabixFormat(TabixFormat.ZERO_BASED, 1, 2, 0, '#', 0);
    }

    @Override
    public BafEvidence decode(final String line) {
        final List<String> tokens = splitter.splitToList(line);
        if (tokens.size() != 4) {
            throw new IllegalArgumentException("Invalid number of columns: " + tokens.size());
        }
        final String contig = tokens.get(0);
        final int position = Integer.parseUnsignedInt(tokens.get(1)) + 1; // Adjust for 0-based indexing
        final double value = Double.parseDouble(tokens.get(2));
        final String sample = tokens.get(3);
        return new BafEvidence(sample, contig, position, value);
    }

    @Override
    public boolean canDecode(final String path) {
        String toDecode = path.toLowerCase();
        if ( IOUtil.hasBlockCompressedExtension(toDecode) ) {
            toDecode = toDecode.substring(0, toDecode.lastIndexOf('.'));
        }
        return toDecode.endsWith(FORMAT_SUFFIX);
    }

    @Override
    public Object readActualHeader(final LineIterator reader) { return null; }

    @Override
    public FeatureOutputStream<BafEvidence> makeSink( final GATKPath path,
                                                      final SAMSequenceDictionary dict,
                                                      final List<String> sampleNames,
                                                      final int compressionLevel ) {
        return new FeatureOutputStream<>(path, getTabixFormat(), BafEvidenceCodec::encode,
                                            dict, compressionLevel);
    }

    @Override
    public void encode( final BafEvidence ev, final FeatureOutputStream<BafEvidence> os ) {
        os.write(ev);
    }

    public static String encode( final BafEvidence ev ) {
        final List<String> columns = Arrays.asList(
                ev.getContig(),
                Integer.toString(ev.getStart() - 1),
                Double.toString(ev.getValue()),
                ev.getSample()
        );
        return String.join(COL_DELIMITER, columns);
    }
}
