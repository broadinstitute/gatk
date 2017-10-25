package org.broadinstitute.hellbender.tools.copynumber.multidimensional.model;

import org.broadinstitute.hellbender.tools.copynumber.formats.collections.SampleLocatableCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SampleMetadata;
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
public final class ModeledSegmentCollection extends SampleLocatableCollection<ModeledSegment> {
    private static final String DOUBLE_FORMAT = "%6.6f";    //TODO replace this with MultidimensionalModeller.DOUBLE_FORMAT from sl_wgs_acnv branch

    enum ModeledSegmentTableColumn {
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
        MINOR_ALLELE_FRACTION_POSTERIOR_90;

        static final TableColumnCollection COLUMNS = new TableColumnCollection((Object[]) values());
    }

    private static final Function<DataLine, ModeledSegment> MODELED_SEGMENT_RECORD_FROM_DATA_LINE_DECODER = dataLine -> {
        final String contig = dataLine.get(ModeledSegmentTableColumn.CONTIG);
        final int start = dataLine.getInt(ModeledSegmentTableColumn.START);
        final int end = dataLine.getInt(ModeledSegmentTableColumn.END);
        final int numPointsCopyRatio = dataLine.getInt(ModeledSegmentTableColumn.NUM_POINTS_COPY_RATIO);
        final int numPointsAlleleFraction = dataLine.getInt(ModeledSegmentTableColumn.NUM_POINTS_ALLELE_FRACTION);
        final double log2CopyRatioPosterior10 = dataLine.getDouble(ModeledSegmentTableColumn.LOG2_COPY_RATIO_POSTERIOR_10);
        final double log2CopyRatioPosterior50 = dataLine.getDouble(ModeledSegmentTableColumn.LOG2_COPY_RATIO_POSTERIOR_50);
        final double log2CopyRatioPosterior90 = dataLine.getDouble(ModeledSegmentTableColumn.LOG2_COPY_RATIO_POSTERIOR_90);
        final double minorAlleleFractionPosterior10 = dataLine.getDouble(ModeledSegmentTableColumn.MINOR_ALLELE_FRACTION_POSTERIOR_10);
        final double minorAlleleFractionPosterior50 = dataLine.getDouble(ModeledSegmentTableColumn.MINOR_ALLELE_FRACTION_POSTERIOR_50);
        final double minorAlleleFractionPosterior90 = dataLine.getDouble(ModeledSegmentTableColumn.MINOR_ALLELE_FRACTION_POSTERIOR_90);
        final SimpleInterval interval = new SimpleInterval(contig, start, end);
        return new ModeledSegment(interval, numPointsCopyRatio, numPointsAlleleFraction, 
                new ModeledSegment.SimplePosteriorSummary(log2CopyRatioPosterior10, log2CopyRatioPosterior50, log2CopyRatioPosterior90),
                new ModeledSegment.SimplePosteriorSummary(minorAlleleFractionPosterior10, minorAlleleFractionPosterior50, minorAlleleFractionPosterior90));
    };

    private static final BiConsumer<ModeledSegment, DataLine> MODELED_SEGMENT_RECORD_TO_DATA_LINE_ENCODER = (modeledSegment, dataLine) ->
            dataLine.append(modeledSegment.getContig())
                    .append(modeledSegment.getStart())
                    .append(modeledSegment.getEnd())
                    .append(modeledSegment.getNumPointsCopyRatio())
                    .append(modeledSegment.getNumPointsAlleleFraction())
                    .append(String.format(DOUBLE_FORMAT, modeledSegment.getLog2CopyRatioSimplePosteriorSummary().getDecile10()))
                    .append(String.format(DOUBLE_FORMAT, modeledSegment.getLog2CopyRatioSimplePosteriorSummary().getDecile50()))
                    .append(String.format(DOUBLE_FORMAT, modeledSegment.getLog2CopyRatioSimplePosteriorSummary().getDecile90()))
                    .append(String.format(DOUBLE_FORMAT, modeledSegment.getMinorAlleleFractionSimplePosteriorSummary().getDecile10()))
                    .append(String.format(DOUBLE_FORMAT, modeledSegment.getMinorAlleleFractionSimplePosteriorSummary().getDecile50()))
                    .append(String.format(DOUBLE_FORMAT, modeledSegment.getMinorAlleleFractionSimplePosteriorSummary().getDecile90()));

    public ModeledSegmentCollection(final File inputFile) {
        super(inputFile, ModeledSegmentTableColumn.COLUMNS, MODELED_SEGMENT_RECORD_FROM_DATA_LINE_DECODER, MODELED_SEGMENT_RECORD_TO_DATA_LINE_ENCODER);
    }

    public ModeledSegmentCollection(final SampleMetadata sampleMetadata,
                                    final List<ModeledSegment> modeledSegments) {
        super(sampleMetadata, modeledSegments, ModeledSegmentTableColumn.COLUMNS, MODELED_SEGMENT_RECORD_FROM_DATA_LINE_DECODER, MODELED_SEGMENT_RECORD_TO_DATA_LINE_ENCODER);
    }
}