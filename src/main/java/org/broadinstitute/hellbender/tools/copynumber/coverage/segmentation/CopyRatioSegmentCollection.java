package org.broadinstitute.hellbender.tools.copynumber.coverage.segmentation;

import org.broadinstitute.hellbender.tools.copynumber.formats.collections.SampleLocatableCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SampleMetadata;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;

import java.io.File;
import java.util.List;
import java.util.function.BiConsumer;
import java.util.function.Function;

public final class CopyRatioSegmentCollection extends SampleLocatableCollection<CopyRatioSegment> {
    enum CopyRatioSegmentTableColumn {
        CONTIG,
        START,
        END,
        NUM_POINTS_COPY_RATIO,
        MEAN_LOG2_COPY_RATIO;

        static final TableColumnCollection COLUMNS = new TableColumnCollection((Object[]) values());
    }
    private static final Function<DataLine, CopyRatioSegment> COPY_RATIO_SEGMENT_RECORD_FROM_DATA_LINE_DECODER = dataLine -> {
        final String contig = dataLine.get(CopyRatioSegmentTableColumn.CONTIG);
        final int start = dataLine.getInt(CopyRatioSegmentTableColumn.START);
        final int end = dataLine.getInt(CopyRatioSegmentTableColumn.END);
        final int numPoints = dataLine.getInt(CopyRatioSegmentTableColumn.NUM_POINTS_COPY_RATIO);
        final double meanLog2CopyRatio = dataLine.getDouble(CopyRatioSegmentTableColumn.MEAN_LOG2_COPY_RATIO);
        final SimpleInterval interval = new SimpleInterval(contig, start, end);
        return new CopyRatioSegment(interval, numPoints, meanLog2CopyRatio);
    };

    private static final BiConsumer<CopyRatioSegment, DataLine> COPY_RATIO_SEGMENT_RECORD_TO_DATA_LINE_ENCODER = (copyRatioSegment, dataLine) ->
            dataLine.append(copyRatioSegment.getInterval().getContig())
                    .append(copyRatioSegment.getInterval().getStart())
                    .append(copyRatioSegment.getInterval().getEnd())
                    .append(copyRatioSegment.getNumPoints())
                    .append(copyRatioSegment.getMeanLog2CopyRatio());

    public CopyRatioSegmentCollection(final File inputFile) {
        super(inputFile, CopyRatioSegmentTableColumn.COLUMNS, COPY_RATIO_SEGMENT_RECORD_FROM_DATA_LINE_DECODER, COPY_RATIO_SEGMENT_RECORD_TO_DATA_LINE_ENCODER);
    }

    public CopyRatioSegmentCollection(final SampleMetadata sampleMetadata,
                                      final List<CopyRatioSegment> copyRatioSegments) {
        super(sampleMetadata, copyRatioSegments, CopyRatioSegmentTableColumn.COLUMNS, COPY_RATIO_SEGMENT_RECORD_FROM_DATA_LINE_DECODER, COPY_RATIO_SEGMENT_RECORD_TO_DATA_LINE_ENCODER);
    }
}