package org.broadinstitute.hellbender.tools.copynumber.formats.collections;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SampleMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.CalledCopyRatioSegment;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.CopyRatioSegment;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;

import java.io.File;
import java.util.Arrays;
import java.util.List;
import java.util.function.BiConsumer;
import java.util.function.Function;

public final class CalledCopyRatioSegmentCollection extends SampleLocatableCollection<CalledCopyRatioSegment> {
    enum CalledCopyRatioSegmentTableColumn {
        CONTIG,
        START,
        END,
        NUM_POINTS_COPY_RATIO,
        MEAN_LOG2_COPY_RATIO,
        CALL;

        static final TableColumnCollection COLUMNS = new TableColumnCollection((Object[]) values());
    }

    private static final Function<DataLine, CalledCopyRatioSegment> CALLED_COPY_RATIO_SEGMENT_RECORD_FROM_DATA_LINE_DECODER = dataLine -> {
        final String contig = dataLine.get(CalledCopyRatioSegmentTableColumn.CONTIG);
        final int start = dataLine.getInt(CalledCopyRatioSegmentTableColumn.START);
        final int end = dataLine.getInt(CalledCopyRatioSegmentTableColumn.END);
        final int numPoints = dataLine.getInt(CalledCopyRatioSegmentTableColumn.NUM_POINTS_COPY_RATIO);
        final double meanLog2CopyRatio = dataLine.getDouble(CalledCopyRatioSegmentTableColumn.MEAN_LOG2_COPY_RATIO);
        final String callOutputString = dataLine.get(CalledCopyRatioSegmentTableColumn.CALL);
        final CalledCopyRatioSegment.Call call = Arrays.stream(CalledCopyRatioSegment.Call.values())
                .filter(c -> c.getOutputString().equals(callOutputString)).findFirst().orElse(null);
        if (call == null) {
            throw new UserException.BadInput(String.format("Invalid call: %s", callOutputString));
        }
        final SimpleInterval interval = new SimpleInterval(contig, start, end);
        return new CalledCopyRatioSegment(new CopyRatioSegment(interval, numPoints, meanLog2CopyRatio), call);
    };

    private static final BiConsumer<CalledCopyRatioSegment, DataLine> CALLED_COPY_RATIO_SEGMENT_RECORD_TO_DATA_LINE_ENCODER = (calledCopyRatioSegment, dataLine) ->
            dataLine.append(calledCopyRatioSegment.getInterval().getContig())
                    .append(calledCopyRatioSegment.getInterval().getStart())
                    .append(calledCopyRatioSegment.getInterval().getEnd())
                    .append(calledCopyRatioSegment.getNumPoints())
                    .append(calledCopyRatioSegment.getMeanLog2CopyRatio())
                    .append(calledCopyRatioSegment.getCall().getOutputString());

    public CalledCopyRatioSegmentCollection(final File inputFile) {
        super(inputFile, CalledCopyRatioSegmentTableColumn.COLUMNS, CALLED_COPY_RATIO_SEGMENT_RECORD_FROM_DATA_LINE_DECODER, CALLED_COPY_RATIO_SEGMENT_RECORD_TO_DATA_LINE_ENCODER);
    }

    public CalledCopyRatioSegmentCollection(final SampleMetadata sampleMetadata,
                                            final List<CalledCopyRatioSegment> calledCopyRatioSegments) {
        super(sampleMetadata, calledCopyRatioSegments, CalledCopyRatioSegmentTableColumn.COLUMNS, CALLED_COPY_RATIO_SEGMENT_RECORD_FROM_DATA_LINE_DECODER, CALLED_COPY_RATIO_SEGMENT_RECORD_TO_DATA_LINE_ENCODER);
    }
}