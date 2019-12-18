package org.broadinstitute.hellbender.utils.codecs;

import com.google.common.base.Splitter;
import htsjdk.samtools.util.IOUtil;
import htsjdk.tribble.AsciiFeatureCodec;
import htsjdk.tribble.readers.LineIterator;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.sv.DepthEvidence;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class DepthEvidenceCodec extends AsciiFeatureCodec<DepthEvidence> {

    public static final String FORMAT_SUFFIX = ".rd.txt";
    public static final String COL_DELIMITER = "\t";
    public static final String HEADER_CHAR = "#";
    private static final Splitter splitter = Splitter.on(COL_DELIMITER);

    public DepthEvidenceCodec() {
        super(DepthEvidence.class);
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
        return new DepthEvidenceMetadata(reader.next());
    }

    public final class DepthEvidenceMetadata {
        final List<String> samples;
        public DepthEvidenceMetadata(final String header) {
            Utils.nonNull(header);
            Utils.validateArg(header.startsWith(HEADER_CHAR), "Expected header starting with " + HEADER_CHAR);
            final String[] tokens = header.split(COL_DELIMITER);
            samples = IntStream.range(3, tokens.length).mapToObj(i -> tokens[i]).collect(Collectors.toList());
        }

        public List<String> getSamples() {
            return samples;
        }
    }
}
