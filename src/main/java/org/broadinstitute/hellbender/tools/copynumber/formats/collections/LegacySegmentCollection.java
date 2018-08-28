package org.broadinstitute.hellbender.tools.copynumber.formats.collections;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SampleLocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.LegacySegment;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;
import java.util.function.BiConsumer;
import java.util.function.Function;
import java.util.stream.Stream;

/**
 * Represents a CBS-style segmentation to enable IGV-compatible plotting.
 *
 * IGV ignores column headers and requires that no other headers are present.
 * We use the conventional CBS-style column headers (which includes the sample name)
 * and suppress the SAM-style metadata header (which breaks the contract for construction from input files).
 * See <a href="http://software.broadinstitute.org/cancer/software/genepattern/file-formats-guide#CBS">
 *     http://software.broadinstitute.org/cancer/software/genepattern/file-formats-guide#CBS</a>
 * and <a href="https://software.broadinstitute.org/software/igv/SEG">
 *     https://software.broadinstitute.org/software/igv/SEG</a>.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class LegacySegmentCollection extends AbstractSampleLocatableCollection<LegacySegment> {
    //note to developers: repeat the column headers in Javadoc so that they are viewable when linked
    /**
     * Sample, Chromosome, Start, End, Num_Probes, Segment_Mean
     */
    enum LegacySegmentTableColumn {
        SAMPLE("Sample"),
        CHROMOSOME("Chromosome"),
        START("Start"),
        END("End"),
        NUM_PROBES("Num_Probes"),
        SEGMENT_MEAN("Segment_Mean");

        private final String columnName;

        LegacySegmentTableColumn(final String columnName) {
            this.columnName = columnName;
        }

        static final TableColumnCollection COLUMNS = new TableColumnCollection(
                Stream.of(values()).map(c -> c.columnName).toArray());
    }

    private static final Function<DataLine, LegacySegment> LEGACY_SEGMENT_DATA_LINE_TO_RECORD_FUNCTION = dataLine -> {
        final String sampleName = dataLine.get(LegacySegmentTableColumn.SAMPLE.columnName);
        final String contig = dataLine.get(LegacySegmentTableColumn.CHROMOSOME.columnName);
        final int start = dataLine.getInt(LegacySegmentTableColumn.START.columnName);
        final int end = dataLine.getInt(LegacySegmentTableColumn.END.columnName);
        final int numProbes = dataLine.getInt(LegacySegmentTableColumn.NUM_PROBES.columnName);
        final double segmentMean = dataLine.getDouble(LegacySegmentTableColumn.SEGMENT_MEAN.columnName);
        final SimpleInterval interval = new SimpleInterval(contig, start, end);
        return new LegacySegment(sampleName, interval, numProbes, segmentMean);
    };

    private static final BiConsumer<LegacySegment, DataLine> LEGACY_SEGMENT_RECORD_AND_DATA_LINE_BI_CONSUMER = (LegacySegment, dataLine) ->
            dataLine.append(LegacySegment.getSampleName())
                    .append(LegacySegment.getContig())
                    .append(LegacySegment.getStart())
                    .append(LegacySegment.getEnd())
                    .append(LegacySegment.getNumProbes())
                    .append(formatDouble(LegacySegment.getSegmentMean()));

    public LegacySegmentCollection(final SampleLocatableMetadata metadata,
                                   final List<LegacySegment> legacySegments) {
        super(metadata, legacySegments, LegacySegmentTableColumn.COLUMNS, LEGACY_SEGMENT_DATA_LINE_TO_RECORD_FUNCTION, LEGACY_SEGMENT_RECORD_AND_DATA_LINE_BI_CONSUMER);
    }

    // output of SAM-style header is suppressed
    @Override
    public void write(final File outputFile) {
        try (final RecordWriter recordWriter = new RecordWriter(new FileWriter(outputFile, true))) {
            recordWriter.writeAllRecords(getRecords());
        } catch (final IOException e) {
            throw new UserException.CouldNotCreateOutputFile(outputFile, e);
        }
    }
}