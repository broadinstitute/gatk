package org.broadinstitute.hellbender.tools.copynumber.formats.collections;

import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SampleLocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.IntegerCopyNumberSegment;
import org.broadinstitute.hellbender.tools.copynumber.gcnv.IntegerCopyNumberState;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;

import java.io.File;
import java.util.List;
import java.util.function.BiConsumer;
import java.util.function.Function;

/**
 * Represents a collection of {@link IntegerCopyNumberSegment} for a sample.
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public final class IntegerCopyNumberSegmentCollection extends AbstractSampleLocatableCollection<IntegerCopyNumberSegment> {
    //note to developers: repeat the column headers in Javadoc so that they are viewable when linked
    /**
     * CONTIG, START, END, NUM_POINTS, CALL_COPY_NUMBER, BASELINE_COPY_NUMBER,
     * QUALITY_SOME_CALLED, QUALITY_ALL_CALLED, QUALITY_START, QUALITY_END
     */
    enum IntegerCopyNumberSegmentTableColumn {
        CONTIG,
        START,
        END,
        NUM_POINTS,
        CALL_COPY_NUMBER,
        BASELINE_COPY_NUMBER,
        QUALITY_SOME_CALLED,
        QUALITY_ALL_CALLED,
        QUALITY_START,
        QUALITY_END;

        static final TableColumnCollection COLUMNS = new TableColumnCollection((Object[]) values());
    }

    private static final Function<DataLine, IntegerCopyNumberSegment> INTEGER_COPY_NUMBER_SEGMENT_RECORD_DECODER = dataLine -> {
        final String contig = dataLine.get(IntegerCopyNumberSegmentTableColumn.CONTIG);
        final int start = dataLine.getInt(IntegerCopyNumberSegmentTableColumn.START);
        final int end = dataLine.getInt(IntegerCopyNumberSegmentTableColumn.END);
        final int numPoints = dataLine.getInt(IntegerCopyNumberSegmentTableColumn.NUM_POINTS);
        final int callCopyNumber = dataLine.getInt(IntegerCopyNumberSegmentTableColumn.CALL_COPY_NUMBER);
        final int baselineCopyNumber = dataLine.getInt(IntegerCopyNumberSegmentTableColumn.BASELINE_COPY_NUMBER);
        final double qualitySomeCalled = dataLine.getDouble(IntegerCopyNumberSegmentTableColumn.QUALITY_SOME_CALLED);
        final double qualityAllCalled = dataLine.getDouble(IntegerCopyNumberSegmentTableColumn.QUALITY_ALL_CALLED);
        final double qualityStart = dataLine.getDouble(IntegerCopyNumberSegmentTableColumn.QUALITY_START);
        final double qualityEnd = dataLine.getDouble(IntegerCopyNumberSegmentTableColumn.QUALITY_END);
        return new IntegerCopyNumberSegment(
                new SimpleInterval(contig, start, end),
                new IntegerCopyNumberState(callCopyNumber),
                new IntegerCopyNumberState(baselineCopyNumber),
                numPoints, qualitySomeCalled, qualityAllCalled, qualityStart, qualityEnd);
    };

    private static final BiConsumer<IntegerCopyNumberSegment, DataLine> INTEGER_COPY_NUMBER_SEGMENT_RECORD_ENCODER =
            (integerCopyNumberSegment, dataLine) ->
                    dataLine.append(integerCopyNumberSegment.getContig())
                            .append(integerCopyNumberSegment.getStart())
                            .append(integerCopyNumberSegment.getEnd())
                            .append(integerCopyNumberSegment.getNumPoints())
                            .append(integerCopyNumberSegment.getCallIntegerCopyNumberState().getCopyNumber())
                            .append(integerCopyNumberSegment.getBaselineIntegerCopyNumberState().getCopyNumber())
                            .append(formatDouble(integerCopyNumberSegment.getQualitySomeCalled()))
                            .append(formatDouble(integerCopyNumberSegment.getQualityAllCalled()))
                            .append(formatDouble(integerCopyNumberSegment.getQualityStart()))
                            .append(formatDouble(integerCopyNumberSegment.getQualityEnd()));

    public IntegerCopyNumberSegmentCollection(final File inputFile) {
        super(inputFile, IntegerCopyNumberSegmentTableColumn.COLUMNS,
                INTEGER_COPY_NUMBER_SEGMENT_RECORD_DECODER, INTEGER_COPY_NUMBER_SEGMENT_RECORD_ENCODER);
    }

    public IntegerCopyNumberSegmentCollection(final SampleLocatableMetadata metadata,
                                              final List<IntegerCopyNumberSegment> integerCopyNumberSegmentList) {
        super(metadata, integerCopyNumberSegmentList, IntegerCopyNumberSegmentTableColumn.COLUMNS,
                INTEGER_COPY_NUMBER_SEGMENT_RECORD_DECODER, INTEGER_COPY_NUMBER_SEGMENT_RECORD_ENCODER);
    }
}
