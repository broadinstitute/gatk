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
     * CONTIG, START, END, NUM_SPANNING_INTERVALS, CALL_COPY_NUMBER, BASELINE_COPY_NUMBER,
     * SOME_QUALITY, EXACT_QUALITY, START_QUALITY, END_QUALITY
     */
    enum IntegerCopyNumberSegmentTableColumn {
        CONTIG,
        START,
        END,
        NUM_SPANNING_INTERVALS,
        CALL_COPY_NUMBER,
        BASELINE_COPY_NUMBER,
        SOME_QUALITY,
        EXACT_QUALITY,
        START_QUALITY,
        END_QUALITY;

        static final TableColumnCollection COLUMNS = new TableColumnCollection((Object[]) values());
    }

    private static final Function<DataLine, IntegerCopyNumberSegment> INTEGER_COPY_NUMBER_SEGMENT_RECORD_DECODER = dataLine -> {
        final String contig = dataLine.get(IntegerCopyNumberSegmentTableColumn.CONTIG);
        final int start = dataLine.getInt(IntegerCopyNumberSegmentTableColumn.START);
        final int end = dataLine.getInt(IntegerCopyNumberSegmentTableColumn.END);
        final int numSpanningIntervals = dataLine.getInt(IntegerCopyNumberSegmentTableColumn.NUM_SPANNING_INTERVALS);
        final int callCopyNumber = dataLine.getInt(IntegerCopyNumberSegmentTableColumn.CALL_COPY_NUMBER);
        final int baselineCopyNumber = dataLine.getInt(IntegerCopyNumberSegmentTableColumn.BASELINE_COPY_NUMBER);
        final double someQuality = dataLine.getDouble(IntegerCopyNumberSegmentTableColumn.SOME_QUALITY);
        final double exactQuality = dataLine.getDouble(IntegerCopyNumberSegmentTableColumn.EXACT_QUALITY);
        final double startQuality = dataLine.getDouble(IntegerCopyNumberSegmentTableColumn.START_QUALITY);
        final double endQuality = dataLine.getDouble(IntegerCopyNumberSegmentTableColumn.END_QUALITY);
        return new IntegerCopyNumberSegment(
                new SimpleInterval(contig, start, end),
                new IntegerCopyNumberState(callCopyNumber),
                new IntegerCopyNumberState(baselineCopyNumber),
                numSpanningIntervals, someQuality, exactQuality, startQuality, endQuality);
    };

    private static final BiConsumer<IntegerCopyNumberSegment, DataLine> INTEGER_COPY_NUMBER_SEGMENT_RECORD_ENCODER =
            (integerCopyNumberSegment, dataLine) ->
                    dataLine.append(integerCopyNumberSegment.getContig())
                            .append(integerCopyNumberSegment.getStart())
                            .append(integerCopyNumberSegment.getEnd())
                            .append(integerCopyNumberSegment.getNumSpanningIntervals())
                            .append(integerCopyNumberSegment.getCallIntegerCopyNumberState().getCopyNumber())
                            .append(integerCopyNumberSegment.getBaselineIntegerCopyNumberState().getCopyNumber())
                            .append(formatDouble(integerCopyNumberSegment.getSomeQuality()))
                            .append(formatDouble(integerCopyNumberSegment.getExactQuality()))
                            .append(formatDouble(integerCopyNumberSegment.getStartQuality()))
                            .append(formatDouble(integerCopyNumberSegment.getEndQuality()));

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