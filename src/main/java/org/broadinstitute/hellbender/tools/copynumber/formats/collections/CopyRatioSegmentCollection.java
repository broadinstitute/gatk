package org.broadinstitute.hellbender.tools.copynumber.formats.collections;

import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SampleLocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.CopyRatioSegment;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;

import java.io.File;
import java.util.List;
import java.util.function.BiConsumer;
import java.util.function.Function;

/**
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class CopyRatioSegmentCollection extends AbstractSampleLocatableCollection<CopyRatioSegment> {
    //note to developers: repeat the column headers in Javadoc so that they are viewable when linked
    /**
     * CONTIG, START, END, NUM_POINTS_COPY_RATIO, MEAN_LOG2_COPY_RATIO
     */
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
                    .append(formatDouble(copyRatioSegment.getMeanLog2CopyRatio()));

    public CopyRatioSegmentCollection(final File inputFile) {
        super(inputFile, CopyRatioSegmentTableColumn.COLUMNS, COPY_RATIO_SEGMENT_RECORD_FROM_DATA_LINE_DECODER, COPY_RATIO_SEGMENT_RECORD_TO_DATA_LINE_ENCODER);
    }

    public CopyRatioSegmentCollection(final SampleLocatableMetadata metadata,
                                      final List<CopyRatioSegment> copyRatioSegments) {
        super(metadata, copyRatioSegments, CopyRatioSegmentTableColumn.COLUMNS, COPY_RATIO_SEGMENT_RECORD_FROM_DATA_LINE_DECODER, COPY_RATIO_SEGMENT_RECORD_TO_DATA_LINE_ENCODER);
    }
}