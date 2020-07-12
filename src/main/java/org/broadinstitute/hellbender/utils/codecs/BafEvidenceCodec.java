package org.broadinstitute.hellbender.utils.codecs;

import com.google.common.base.Splitter;
import htsjdk.tribble.AsciiFeatureCodec;
import htsjdk.tribble.readers.LineIterator;
import org.broadinstitute.hellbender.tools.sv.BafEvidence;

import java.util.List;

public class BafEvidenceCodec extends AsciiFeatureCodec<BafEvidence> {

    public static final String FORMAT_SUFFIX = ".BAF.txt.gz";
    public static final String COL_DELIMITER = "\t";
    private static final Splitter splitter = Splitter.on(COL_DELIMITER);

    public BafEvidenceCodec() {
        super(BafEvidence.class);
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
        return path.endsWith(FORMAT_SUFFIX);
    }

    @Override
    public Object readActualHeader(final LineIterator reader) { return null; }
}
