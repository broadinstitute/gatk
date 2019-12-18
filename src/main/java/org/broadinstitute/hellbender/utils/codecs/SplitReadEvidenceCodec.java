package org.broadinstitute.hellbender.utils.codecs;

import com.google.common.base.Splitter;
import htsjdk.tribble.AsciiFeatureCodec;
import htsjdk.tribble.index.tabix.TabixFormat;
import htsjdk.tribble.readers.LineIterator;
import org.broadinstitute.hellbender.tools.sv.SplitReadEvidence;

import java.util.List;

public class SplitReadEvidenceCodec extends AsciiFeatureCodec<SplitReadEvidence> {

    public static final String FORMAT_SUFFIX = ".txt.gz";
    public static final char COL_DELIMITER = '\t';
    public static final String DIRECTION_RIGHT = "right";
    public static final String DIRECTION_LEFT = "left";
    private static final Splitter splitter = Splitter.on(COL_DELIMITER);

    public SplitReadEvidenceCodec() {
        super(SplitReadEvidence.class);
    }

    @Override
    public SplitReadEvidence decode(final String line) {
        final List<String> tokens = splitter.splitToList(line);
        if (tokens.size() != 5) {
            throw new IllegalArgumentException("Invalid number of columns: " + tokens.size());
        }
        final String contig = tokens.get(0);
        final int position = Integer.parseUnsignedInt(tokens.get(1)) + 1; // Adjust for 0-based indexing
        if (!tokens.get(2).equals(DIRECTION_LEFT) && !tokens.get(2).equals(DIRECTION_RIGHT)) {
            throw new IllegalArgumentException("Invalid direction: " + tokens.get(2));
        }
        final boolean strand = tokens.get(2).equals(DIRECTION_RIGHT);
        final int count = Integer.parseUnsignedInt(tokens.get(3));
        final String sample = tokens.get(4);
        return new SplitReadEvidence(sample, contig, position, count, strand);
    }


    @Override
    public TabixFormat getTabixFormat() {
        return TabixFormat.BED;
    }

    @Override
    public boolean canDecode(final String path) {
        return path.endsWith(FORMAT_SUFFIX);
    } // TODO

    @Override
    public Object readActualHeader(final LineIterator reader) { return null; }
}
