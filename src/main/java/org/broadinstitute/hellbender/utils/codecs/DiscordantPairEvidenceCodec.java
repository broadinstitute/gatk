package org.broadinstitute.hellbender.utils.codecs;

import com.google.common.base.Splitter;
import htsjdk.samtools.util.IOUtil;
import htsjdk.tribble.AsciiFeatureCodec;
import htsjdk.tribble.readers.LineIterator;
import org.broadinstitute.hellbender.tools.sv.DiscordantPairEvidence;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;

import java.util.List;

public class DiscordantPairEvidenceCodec extends AsciiFeatureCodec<DiscordantPairEvidence> {

    public static final String FORMAT_SUFFIX = ".pe.txt";
    public static final String COL_DELIMITER = "\t";
    private static final Splitter splitter = Splitter.on(COL_DELIMITER);

    public DiscordantPairEvidenceCodec() {
        super(DiscordantPairEvidence.class);
    }

    @Override
    public DiscordantPairEvidence decode(final String line) {
        final List<String> tokens = splitter.splitToList(line);
        if (tokens.size() != 7) {
            throw new IllegalArgumentException("Invalid number of columns: " + tokens.size());
        }
        final String startContig = tokens.get(0);
        final int start = Integer.parseUnsignedInt(tokens.get(1)) + 1; // Adjust for 0-based indexing
        final boolean startStrand = tokens.get(2).equals(SVCallRecord.STRAND_PLUS);
        final String endContig = tokens.get(3);
        final int end = Integer.parseUnsignedInt(tokens.get(4)) + 1; // Adjust for 0-based indexing
        final boolean endStrand = tokens.get(5).equals(SVCallRecord.STRAND_PLUS);
        final String sample = tokens.get(6);
        return new DiscordantPairEvidence(sample, startContig, start, startStrand, endContig, end, endStrand);
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
    public Object readActualHeader(final LineIterator reader) { return null; }
}
