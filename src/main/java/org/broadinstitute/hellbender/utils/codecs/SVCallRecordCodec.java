package org.broadinstitute.hellbender.utils.codecs;

import com.google.common.base.Splitter;
import htsjdk.tribble.AsciiFeatureCodec;
import htsjdk.tribble.index.tabix.TabixFormat;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.StructuralVariantType;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class SVCallRecordCodec extends AsciiFeatureCodec<SVCallRecord> {

    public static final String FORMAT_SUFFIX = ".tsv.gz";
    public static final String COL_DELIMITER = "\t";
    public static String STRAND_PLUS = "+";
    public static String STRAND_MINUS = "-";
    private static final Splitter splitter = Splitter.on(COL_DELIMITER);

    public SVCallRecordCodec() {
        super(SVCallRecord.class);
    }

    @Override
    public SVCallRecord decode(final String line) {
        final List<String> tokens = splitter.splitToList(line);
        if (tokens.size() != 10) {
            throw new IllegalArgumentException("Invalid number of columns: " + tokens.size());
        }
        return new SVCallRecord(
                tokens.get(0),
                Integer.parseUnsignedInt(tokens.get(1)) + 1, // Convert to 1-based indexing
                tokens.get(2).equals(STRAND_PLUS),
                tokens.get(3),
                Integer.parseUnsignedInt(tokens.get(4)) + 1,
                tokens.get(5).equals(STRAND_PLUS),
                StructuralVariantType.valueOf(tokens.get(6)),
                Integer.parseInt(tokens.get(7)),
                Arrays.asList(tokens.get(8)),
                Collections.singleton(tokens.get(9))
        );
    }

    @Override
    public TabixFormat getTabixFormat() {
        return TabixFormat.BED;
    }

    @Override
    public boolean canDecode(final String path) {
        return path.endsWith(FORMAT_SUFFIX);
    }

    @Override
    public Object readActualHeader(final LineIterator reader) { return null; }

    public String encode(final SVCallRecord record) {
        final List<String> data = Arrays.asList(
                record.getContig(),
                Integer.toString(record.getStart() - 1), // Convert to 0-based indexing
                record.getStartStrand() ? STRAND_PLUS : STRAND_MINUS,
                record.getEndContig(),
                Integer.toString(record.getEnd() - 1),
                record.getEndStrand() ? STRAND_PLUS : STRAND_MINUS,
                record.getType().name(),
                Integer.toString(record.getLength()),
                String.join(",", record.getAlgorithms()),
                String.join(",", record.getSamples())
        );
        return String.join(COL_DELIMITER, data);
    }
}
