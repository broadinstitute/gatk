package org.broadinstitute.hellbender.utils.codecs;

import com.google.common.base.Splitter;
import htsjdk.samtools.util.IOUtil;
import htsjdk.tribble.AsciiFeatureCodec;
import htsjdk.tribble.index.tabix.TabixFormat;
import htsjdk.tribble.readers.LineIterator;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.sv.DepthEvidence;

import java.util.ArrayList;
import java.util.List;

public class DepthEvidenceCodec extends AsciiFeatureCodec<DepthEvidence> {

    public static final String FORMAT_SUFFIX = ".rd.txt";
    public static final String COL_DELIMITER = "\t";
    private static final Splitter splitter = Splitter.on(COL_DELIMITER);

    public DepthEvidenceCodec() {
        super(DepthEvidence.class);
    }

    @Override
    public TabixFormat getTabixFormat() {
        return new TabixFormat(TabixFormat.ZERO_BASED, 1, 2, 3, '#', 1);
    }

    @Override
    public DepthEvidence decode(final String line) {
        final List<String> tokens = splitter.splitToList(line);
        if (tokens.size() < 3) {
            throw new IllegalArgumentException("Expected at least 3 columns but found " + tokens.size());
        }
        final String contig = tokens.get(0);
        final int start = Integer.parseUnsignedInt(tokens.get(1)) + 1; // Adjust for 0-based indexing
        final int end = Integer.parseUnsignedInt(tokens.get(2)); // Adjust for 0-based indexing and inclusive intervals
        final int numCounts = tokens.size() - 3;
        final int[] counts = new int[numCounts];
        for (int i = 3; i < tokens.size(); i++) {
            counts[i - 3] = Integer.parseUnsignedInt(tokens.get(i));
        }
        return new DepthEvidence(contig, start, end, counts);
    }

    @Override
    public boolean canDecode(final String path) {
        final String toDecode;
        if (IOUtil.hasBlockCompressedExtension(path)) {
            toDecode = path.substring(0, path.lastIndexOf("."));
        } else {
            toDecode = path;
        }
        return toDecode.toLowerCase().endsWith(FORMAT_SUFFIX);
    }

    @Override
    public Object readActualHeader(final LineIterator reader) {
        if (!reader.hasNext()) {
            throw new UserException.BadInput("Depth evidence file did not have a header line");
        }
        return reader.next();
    }

    public static String encode(final DepthEvidence ev) {
        final int[] counts = ev.getCounts();
        final int numCounts = counts.length;
        final List<String> columns = new ArrayList<>(3 + numCounts);
        columns.add(ev.getContig());
        columns.add(Integer.toString(ev.getStart() - 1));
        columns.add(Integer.toString(ev.getEnd()));
        for (int i = 0; i < numCounts; i++) {
            columns.add(Integer.toString(counts[i]));
        }
        return String.join(COL_DELIMITER, columns);
    }
}
