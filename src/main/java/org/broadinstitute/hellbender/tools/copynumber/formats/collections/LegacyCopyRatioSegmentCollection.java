package org.broadinstitute.hellbender.tools.copynumber.formats.collections;

import org.broadinstitute.hellbender.tools.copynumber.formats.records.LegacyCopyRatioSegment;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;

import java.io.File;
import java.util.function.BiConsumer;
import java.util.function.Function;

public class LegacyCopyRatioSegmentCollection extends AbstractSampleLocatableCollection<LegacyCopyRatioSegment> {
    /**
     * SAMPLE, CONTIG, START, END, NUM_POINTS_COPY_RATIO, MEAN_LOG2_COPY_RATIO
     */
    enum LegacyCopyRatioSegmentTableColumn {
        SAMPLE,
        CONTIG,
        START,
        END,
        NUM_POINTS_COPY_RATIO,
        MEAN_LOG2_COPY_RATIO;

        static final TableColumnCollection COLUMNS = new TableColumnCollection((Object[]) values());
    }

    private static final Function<DataLine, LegacyCopyRatioSegment> LEGACY_COPY_RATIO_SEGMENT_RECORD_FROM_DATA_LINE_DECODER = dataLine -> {
        final String sampleName = dataLine.get(LegacyCopyRatioSegmentCollection.LegacyCopyRatioSegmentTableColumn.SAMPLE);
        final String contig = dataLine.get(LegacyCopyRatioSegmentCollection.LegacyCopyRatioSegmentTableColumn.CONTIG);
        final int start = dataLine.getInt(LegacyCopyRatioSegmentCollection.LegacyCopyRatioSegmentTableColumn.START);
        final int end = dataLine.getInt(LegacyCopyRatioSegmentCollection.LegacyCopyRatioSegmentTableColumn.END);
        final int numPoints = dataLine.getInt(LegacyCopyRatioSegmentCollection.LegacyCopyRatioSegmentTableColumn.NUM_POINTS_COPY_RATIO);
        final double meanLog2CopyRatio = dataLine.getDouble(LegacyCopyRatioSegmentCollection.LegacyCopyRatioSegmentTableColumn.MEAN_LOG2_COPY_RATIO);
        final SimpleInterval interval = new SimpleInterval(contig, start, end);
        return new LegacyCopyRatioSegment(sampleName, interval, numPoints, meanLog2CopyRatio);
    };

    private static final BiConsumer<LegacyCopyRatioSegment, DataLine> LEGACY_COPY_RATIO_SEGMENT_RECORD_TO_DATA_LINE_ENCODER = (legacyCopyRatioSegment, dataLine) ->
            dataLine.append(legacyCopyRatioSegment.getSampleName())
                    .append(legacyCopyRatioSegment.getInterval().getContig())
                    .append(legacyCopyRatioSegment.getInterval().getStart())
                    .append(legacyCopyRatioSegment.getInterval().getEnd())
                    .append(legacyCopyRatioSegment.getNumPoints())
                    .append(formatDouble(legacyCopyRatioSegment.getMeanLog2CopyRatio()));

    public LegacyCopyRatioSegmentCollection(final File inputFile) {
        super(inputFile, LegacyCopyRatioSegmentCollection.LegacyCopyRatioSegmentTableColumn.COLUMNS, LEGACY_COPY_RATIO_SEGMENT_RECORD_FROM_DATA_LINE_DECODER, LEGACY_COPY_RATIO_SEGMENT_RECORD_TO_DATA_LINE_ENCODER);
    }

    public LegacyCopyRatioSegmentCollection(final File inputFile, final TableColumnCollection mandatoryColumns, final Function<DataLine, LegacyCopyRatioSegment> recordFromDataLineDecoder, final BiConsumer<LegacyCopyRatioSegment, DataLine> recordToDataLineEncoder) {
        super(inputFile, mandatoryColumns, recordFromDataLineDecoder, recordToDataLineEncoder);
    }
}
