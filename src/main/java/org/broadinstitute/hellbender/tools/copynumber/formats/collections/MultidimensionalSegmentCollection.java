package org.broadinstitute.hellbender.tools.copynumber.formats.collections;

import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SampleLocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.MultidimensionalSegment;
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
public final class MultidimensionalSegmentCollection extends AbstractSampleLocatableCollection<MultidimensionalSegment> {
    //note to developers: repeat the column headers in Javadoc so that they are viewable when linked
    /**
     * CONTIG, START, END, NUM_POINTS_COPY_RATIO, NUM_POINTS_ALLELE_FRACTION, MEAN_LOG2_COPY_RATIO
     */
    enum MultidimensionalSegmentTableColumn {
        CONTIG,
        START,
        END,
        NUM_POINTS_COPY_RATIO,
        NUM_POINTS_ALLELE_FRACTION,
        MEAN_LOG2_COPY_RATIO;

        static final TableColumnCollection COLUMNS = new TableColumnCollection((Object[]) values());
    }

    private static final Function<DataLine, MultidimensionalSegment> MULTIDIMENSIONAL_SEGMENT_DATA_LINE_TO_RECORD_FUNCTION = dataLine -> {
        final String contig = dataLine.get(MultidimensionalSegmentTableColumn.CONTIG);
        final int start = dataLine.getInt(MultidimensionalSegmentTableColumn.START);
        final int end = dataLine.getInt(MultidimensionalSegmentTableColumn.END);
        final int numPointsCopyRatio = dataLine.getInt(MultidimensionalSegmentTableColumn.NUM_POINTS_COPY_RATIO);
        final int numPointsAlleleFraction = dataLine.getInt(MultidimensionalSegmentTableColumn.NUM_POINTS_ALLELE_FRACTION);
        final double meanLog2CopyRatio = dataLine.getDouble(MultidimensionalSegmentTableColumn.MEAN_LOG2_COPY_RATIO);
        final SimpleInterval interval = new SimpleInterval(contig, start, end);
        return new MultidimensionalSegment(interval, numPointsCopyRatio, numPointsAlleleFraction, meanLog2CopyRatio);
    };

    private static final BiConsumer<MultidimensionalSegment, DataLine> MULTIDIMENSIONAL_SEGMENT_RECORD_AND_DATA_LINE_BI_CONSUMER = (alleleFractionSegment, dataLine) ->
            dataLine.append(alleleFractionSegment.getContig())
                    .append(alleleFractionSegment.getStart())
                    .append(alleleFractionSegment.getEnd())
                    .append(alleleFractionSegment.getNumPointsCopyRatio())
                    .append(alleleFractionSegment.getNumPointsAlleleFraction())
                    .append(formatDouble(alleleFractionSegment.getMeanLog2CopyRatio()));

    public MultidimensionalSegmentCollection(final File inputFile) {
        super(inputFile, MultidimensionalSegmentTableColumn.COLUMNS, MULTIDIMENSIONAL_SEGMENT_DATA_LINE_TO_RECORD_FUNCTION, MULTIDIMENSIONAL_SEGMENT_RECORD_AND_DATA_LINE_BI_CONSUMER);
    }

    public MultidimensionalSegmentCollection(final SampleLocatableMetadata metadata,
                                             final List<MultidimensionalSegment> multidimensionalSegments) {
        super(metadata, multidimensionalSegments, MultidimensionalSegmentTableColumn.COLUMNS, MULTIDIMENSIONAL_SEGMENT_DATA_LINE_TO_RECORD_FUNCTION, MULTIDIMENSIONAL_SEGMENT_RECORD_AND_DATA_LINE_BI_CONSUMER);
    }
}