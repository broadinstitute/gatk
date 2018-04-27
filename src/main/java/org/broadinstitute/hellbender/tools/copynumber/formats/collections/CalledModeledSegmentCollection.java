package org.broadinstitute.hellbender.tools.copynumber.formats.collections;

import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SampleLocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.CalledModeledSegment;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;

import java.io.File;
import java.util.List;
import java.util.function.BiConsumer;
import java.util.function.Function;

/**
 * @author Marton Kanasz-Nagy &lt;mkanaszn@broadinstitute.org&gt;
 */
public final class CalledModeledSegmentCollection extends AbstractSampleLocatableCollection<CalledModeledSegment> {
    //note to developers: repeat the column headers in Javadoc so that they are viewable when linked
    /**
     * CONTIG, START, END, NUM_POINTS_COPY_RATIO, NUM_POINTS_ALLELE_FRACTION,
     * LOG2_COPY_RATIO_POSTERIOR_10, LOG2_COPY_RATIO_POSTERIOR_50, LOG2_COPY_RATIO_POSTERIOR_90,
     * MINOR_ALLELE_FRACTION_POSTERIOR_10, MINOR_ALLELE_FRACTION_POSTERIOR_50, MINOR_ALLELE_FRACTION_POSTERIOR_90,
     * CALL_NORMAL, PHRED_SCORE_NORMAL
     */
    enum CalledModeledSegmentTableColumn {
        CONTIG,
        START,
        END,
        NUM_POINTS_COPY_RATIO,
        NUM_POINTS_ALLELE_FRACTION,
        LOG2_COPY_RATIO_POSTERIOR_10,
        LOG2_COPY_RATIO_POSTERIOR_50,
        LOG2_COPY_RATIO_POSTERIOR_90,
        MINOR_ALLELE_FRACTION_POSTERIOR_10,
        MINOR_ALLELE_FRACTION_POSTERIOR_50,
        MINOR_ALLELE_FRACTION_POSTERIOR_90,
        CALL_NORMAL,
        PHRED_SCORE_NORMAL;

        static final TableColumnCollection COLUMNS = new TableColumnCollection((Object[]) values());
    }

    private static final Function<DataLine, CalledModeledSegment> CALLED_MODELED_SEGMENT_RECORD_FROM_DATA_LINE_DECODER = dataLine -> {
        final String contig = dataLine.get(CalledModeledSegmentTableColumn.CONTIG);
        final int start = dataLine.getInt(CalledModeledSegmentTableColumn.START);
        final int end = dataLine.getInt(CalledModeledSegmentTableColumn.END);
        final int numPointsCopyRatio = dataLine.getInt(CalledModeledSegmentTableColumn.NUM_POINTS_COPY_RATIO);
        final int numPointsAlleleFraction = dataLine.getInt(CalledModeledSegmentTableColumn.NUM_POINTS_ALLELE_FRACTION);
        final double log2CopyRatioPosterior10 = dataLine.getDouble(CalledModeledSegmentTableColumn.LOG2_COPY_RATIO_POSTERIOR_10);
        final double log2CopyRatioPosterior50 = dataLine.getDouble(CalledModeledSegmentTableColumn.LOG2_COPY_RATIO_POSTERIOR_50);
        final double log2CopyRatioPosterior90 = dataLine.getDouble(CalledModeledSegmentTableColumn.LOG2_COPY_RATIO_POSTERIOR_90);
        final double minorAlleleFractionPosterior10 = dataLine.getDouble(CalledModeledSegmentTableColumn.MINOR_ALLELE_FRACTION_POSTERIOR_10);
        final double minorAlleleFractionPosterior50 = dataLine.getDouble(CalledModeledSegmentTableColumn.MINOR_ALLELE_FRACTION_POSTERIOR_50);
        final double minorAlleleFractionPosterior90 = dataLine.getDouble(CalledModeledSegmentTableColumn.MINOR_ALLELE_FRACTION_POSTERIOR_90);
        final SimpleInterval interval = new SimpleInterval(contig, start, end);
        final String callNormal = dataLine.get(CalledModeledSegmentTableColumn.CALL_NORMAL);
        final double PHREDScoreNormal = dataLine.getDouble(CalledModeledSegmentTableColumn.PHRED_SCORE_NORMAL);

        return new CalledModeledSegment(interval, numPointsCopyRatio, numPointsAlleleFraction,
                new CalledModeledSegment.SimplePosteriorSummary(log2CopyRatioPosterior10, log2CopyRatioPosterior50, log2CopyRatioPosterior90),
                new CalledModeledSegment.SimplePosteriorSummary(minorAlleleFractionPosterior10, minorAlleleFractionPosterior50, minorAlleleFractionPosterior90),
                callNormal, PHREDScoreNormal);
    };

    private static final BiConsumer<CalledModeledSegment, DataLine> CALLED_MODELED_SEGMENT_RECORD_TO_DATA_LINE_ENCODER = (calledModeledSegment, dataLine) ->
            dataLine.append(calledModeledSegment.getContig())
                    .append(calledModeledSegment.getStart())
                    .append(calledModeledSegment.getEnd())
                    .append(calledModeledSegment.getNumPointsCopyRatio())
                    .append(calledModeledSegment.getNumPointsAlleleFraction())
                    .append(formatDouble(calledModeledSegment.getLog2CopyRatioSimplePosteriorSummary().getDecile10()))
                    .append(formatDouble(calledModeledSegment.getLog2CopyRatioSimplePosteriorSummary().getDecile50()))
                    .append(formatDouble(calledModeledSegment.getLog2CopyRatioSimplePosteriorSummary().getDecile90()))
                    .append(formatDouble(calledModeledSegment.getMinorAlleleFractionSimplePosteriorSummary().getDecile10()))
                    .append(formatDouble(calledModeledSegment.getMinorAlleleFractionSimplePosteriorSummary().getDecile50()))
                    .append(formatDouble(calledModeledSegment.getMinorAlleleFractionSimplePosteriorSummary().getDecile90()))
                    .append(calledModeledSegment.getCallNormal())
                    .append(formatDouble(calledModeledSegment.getPHREDScoreNormal()));

    public CalledModeledSegmentCollection(final File inputFile) {
        super(inputFile, CalledModeledSegmentTableColumn.COLUMNS, CALLED_MODELED_SEGMENT_RECORD_FROM_DATA_LINE_DECODER, CALLED_MODELED_SEGMENT_RECORD_TO_DATA_LINE_ENCODER);
    }

    public CalledModeledSegmentCollection(final SampleLocatableMetadata metadata,
                                    final List<CalledModeledSegment> calledModeledSegments) {
        super(metadata, calledModeledSegments, CalledModeledSegmentTableColumn.COLUMNS, CALLED_MODELED_SEGMENT_RECORD_FROM_DATA_LINE_DECODER, CALLED_MODELED_SEGMENT_RECORD_TO_DATA_LINE_ENCODER);
    }
}