package org.broadinstitute.hellbender.tools.copynumber.formats.collections;

import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SampleLocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.CalledCopyRatioSegment;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.CalledLegacySegment;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
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
 */
public final class CalledLegacySegmentCollection extends AbstractSampleLocatableCollection<CalledLegacySegment> {
    //note to developers: repeat the column headers in Javadoc so that they are viewable when linked
    /**
     * Sample, Chromosome, Start, End, Num_Probes, Call, Segment_Mean
     */
    enum CalledLegacySegmentTableColumn {
        SAMPLE("Sample"),
        CHROMOSOME("Chromosome"),
        START("Start"),
        END("End"),
        NUM_PROBES("Num_Probes"),
        CALL("Call"),
        SEGMENT_MEAN("Segment_Mean");

        private final String columnName;

        CalledLegacySegmentTableColumn(final String columnName) {
            this.columnName = columnName;
        }

        static final TableColumnCollection COLUMNS = new TableColumnCollection(
                Stream.of(values()).map(c -> c.columnName).toArray());
    }

    private static final Function<DataLine, CalledLegacySegment> CALLED_LEGACY_SEGMENT_DATA_LINE_TO_RECORD_FUNCTION = dataLine -> {
        final String sampleName = dataLine.get(CalledLegacySegmentTableColumn.SAMPLE.columnName);
        final String contig = dataLine.get(CalledLegacySegmentTableColumn.CHROMOSOME.columnName);
        final int start = dataLine.getInt(CalledLegacySegmentTableColumn.START.columnName);
        final int end = dataLine.getInt(CalledLegacySegmentTableColumn.END.columnName);
        final int numProbes = dataLine.getInt(CalledLegacySegmentTableColumn.NUM_PROBES.columnName);
        final double segmentMean = dataLine.getDouble(CalledLegacySegmentTableColumn.SEGMENT_MEAN.columnName);
        final String callOutputString = dataLine.get(CalledCopyRatioSegmentCollection.CalledCopyRatioSegmentTableColumn.CALL);
        final CalledCopyRatioSegment.Call call = Arrays.stream(CalledCopyRatioSegment.Call.values())
                .filter(c -> c.getOutputString().equals(callOutputString)).findFirst().orElseThrow(
                        () -> new UserException("Attempting to read an invalid value for " +
                                CalledLegacySegmentTableColumn.CALL +": " + callOutputString +
                                ".  Valid values are " + StringUtils.join(CalledCopyRatioSegment.Call.values(), ", ")
                        ));
        final SimpleInterval interval = new SimpleInterval(contig, start, end);
        return new CalledLegacySegment(sampleName, interval, numProbes, segmentMean, call);
    };

    private static final BiConsumer<CalledLegacySegment, DataLine> CALLED_LEGACY_SEGMENT_RECORD_AND_DATA_LINE_BI_CONSUMER = (calledLegacySegment, dataLine) ->
            dataLine.append(calledLegacySegment.getSampleName())
                    .append(calledLegacySegment.getContig())
                    .append(calledLegacySegment.getStart())
                    .append(calledLegacySegment.getEnd())
                    .append(calledLegacySegment.getNumProbes())
                    .append(calledLegacySegment.getCall().getOutputString())
                    .append(formatDouble(calledLegacySegment.getSegmentMean()));

    public CalledLegacySegmentCollection(final SampleLocatableMetadata metadata,
                                         final List<CalledLegacySegment> calledLegacySegments) {
        super(metadata, calledLegacySegments, CalledLegacySegmentTableColumn.COLUMNS, CALLED_LEGACY_SEGMENT_DATA_LINE_TO_RECORD_FUNCTION, CALLED_LEGACY_SEGMENT_RECORD_AND_DATA_LINE_BI_CONSUMER);
    }


    public CalledLegacySegmentCollection(final File inputFile) {
        super(inputFile, CopyRatioSegmentCollection.CopyRatioSegmentTableColumn.COLUMNS, CALLED_LEGACY_SEGMENT_DATA_LINE_TO_RECORD_FUNCTION, CALLED_LEGACY_SEGMENT_RECORD_AND_DATA_LINE_BI_CONSUMER);
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