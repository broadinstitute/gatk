package org.broadinstitute.hellbender.utils.codecs.copynumber;

import com.google.common.collect.ImmutableList;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMTextHeaderCodec;
import htsjdk.samtools.util.BufferedLineReader;
import htsjdk.tribble.AsciiFeatureCodec;
import htsjdk.tribble.index.tabix.TabixFormat;
import htsjdk.tribble.readers.LineIterator;
import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.formats.CopyNumberFormatsUtils;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.SimpleCountCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.Metadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.MetadataUtils;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SampleLocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.SimpleCount;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.tsv.TableUtils;

import java.util.ArrayList;
import java.util.List;


public final class SimpleCountCodec extends AsciiFeatureCodec<SimpleCount> {

    private static final int SAM_HEADER_LINES_INITIAL_CAPACITY = 10_000;

    private static final String COLUMN_HEADER_STRING = String.join(
            TableUtils.COLUMN_SEPARATOR_STRING,
            SimpleCountCollection.SimpleCountTableColumn.COLUMNS.names());

    private static final int TABIX_FORMAT_SEQUENCE_COLUMN = 1;
    private static final int TABIX_FORMAT_START_POSITION_COLUMN = 2;
    private static final int TABIX_FORMAT_END_POSITION_COLUMN = 3;
    private static final char TABIX_FORMAT_META_CHARACTER = CopyNumberFormatsUtils.COMMENT_PREFIX.charAt(0);

    public static final List<String> SIMPLE_COUNT_CODEC_EXTENSIONS = ImmutableList.of(".counts.tsv", ".counts.tsv.gz");

    public SimpleCountCodec() {
        super(SimpleCount.class);
    }

    @Override
    public SimpleCount decode(final String line) {
        if (line.startsWith(CopyNumberFormatsUtils.COMMENT_PREFIX) || line.startsWith(COLUMN_HEADER_STRING)) {
            return null;
        } else {
            final String[] split = line.split(TableUtils.COLUMN_SEPARATOR_STRING);
            try {
                return new SimpleCount(new SimpleInterval(split[0], Integer.parseInt(split[1]), Integer.parseInt(split[2])), Integer.parseInt(split[3]));
            } catch (final NumberFormatException e) {
                throw new UserException.MalformedFile("Line = " + line + " is not formatted correctly.");
            }
        }
    }

    @Override
    public SampleLocatableMetadata readActualHeader(final LineIterator reader) {
        final List<String> samHeaderLines = new ArrayList<>(SAM_HEADER_LINES_INITIAL_CAPACITY);
        //we check that the SAM header lines and the column header line are present in the correct order, then return the mandatory column header
        boolean isSAMHeaderPresent = false;
        while (reader.hasNext()) {
            final String line = reader.peek();
            if (line.startsWith(CopyNumberFormatsUtils.COMMENT_PREFIX)) {
                isSAMHeaderPresent = true;
                samHeaderLines.add(line);
                reader.next();
            } else {
                if (!isSAMHeaderPresent) {
                    throw new UserException.MalformedFile("SAM header lines must be at the beginning of the file.");
                } else if (!line.startsWith(COLUMN_HEADER_STRING)) {
                    throw new UserException.MalformedFile("File does not have a column header.");
                } else {
                    //we just peeked at the column header line, so we need to advance past it
                    reader.next();
                    break;
                }
            }
        }
        final SAMFileHeader samFileHeader = new SAMTextHeaderCodec()
                .decode(BufferedLineReader.fromString(StringUtils.join(samHeaderLines, System.lineSeparator())), null);
        return MetadataUtils.fromHeader(samFileHeader, Metadata.Type.SAMPLE_LOCATABLE);
    }

    @Override
    public boolean canDecode(final String path) {
        return SIMPLE_COUNT_CODEC_EXTENSIONS.stream().anyMatch(path::endsWith);
    }

    @Override
    public TabixFormat getTabixFormat() {
        return new TabixFormat(
                TabixFormat.GENERIC_FLAGS,
                TABIX_FORMAT_SEQUENCE_COLUMN,
                TABIX_FORMAT_START_POSITION_COLUMN,
                TABIX_FORMAT_END_POSITION_COLUMN,
                TABIX_FORMAT_META_CHARACTER,
                0);
    }

    public static String encode( final SimpleCount simpleCount ) {
        return simpleCount.getContig() + "\t" +
                simpleCount.getStart() + "\t" +
                simpleCount.getEnd() + "\t" +
                simpleCount.getCount();
    }
}
