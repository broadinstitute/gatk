package org.broadinstitute.hellbender.tools.copynumber.formats.collections;

import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SampleMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.AlleleFractionSegment;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;

import java.io.File;
import java.util.List;
import java.util.function.BiConsumer;
import java.util.function.Function;

/**
 * Represents an allele-fraction segmentation.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class AlleleFractionSegmentCollection extends SampleLocatableCollection<AlleleFractionSegment> {
    enum AlleleFractionSegmentTableColumn {
        CONTIG,
        START,
        END,
        NUM_POINTS_ALLELE_FRACTION;

        static final TableColumnCollection COLUMNS = new TableColumnCollection((Object[]) values());
    }

    private static final Function<DataLine, AlleleFractionSegment> ALLELE_FRACTION_SEGMENT_DATA_LINE_TO_RECORD_FUNCTION = dataLine -> {
        final String contig = dataLine.get(AlleleFractionSegmentTableColumn.CONTIG);
        final int start = dataLine.getInt(AlleleFractionSegmentTableColumn.START);
        final int end = dataLine.getInt(AlleleFractionSegmentTableColumn.END);
        final int numPoints = dataLine.getInt(AlleleFractionSegmentTableColumn.NUM_POINTS_ALLELE_FRACTION);
        final SimpleInterval interval = new SimpleInterval(contig, start, end);
        return new AlleleFractionSegment(interval, numPoints);
    };

    private static final BiConsumer<AlleleFractionSegment, DataLine> ALLELE_FRACTION_SEGMENT_RECORD_AND_DATA_LINE_BI_CONSUMER = (alleleFractionSegment, dataLine) ->
            dataLine.append(alleleFractionSegment.getContig())
                    .append(alleleFractionSegment.getStart())
                    .append(alleleFractionSegment.getEnd())
                    .append(alleleFractionSegment.getNumPoints());

    public AlleleFractionSegmentCollection(final File inputFile) {
        super(inputFile, AlleleFractionSegmentTableColumn.COLUMNS, ALLELE_FRACTION_SEGMENT_DATA_LINE_TO_RECORD_FUNCTION, ALLELE_FRACTION_SEGMENT_RECORD_AND_DATA_LINE_BI_CONSUMER);
    }

    public AlleleFractionSegmentCollection(final SampleMetadata sampleMetadata,
                                           final List<AlleleFractionSegment> AlleleFractionSegments) {
        super(sampleMetadata, AlleleFractionSegments, AlleleFractionSegmentTableColumn.COLUMNS, ALLELE_FRACTION_SEGMENT_DATA_LINE_TO_RECORD_FUNCTION, ALLELE_FRACTION_SEGMENT_RECORD_AND_DATA_LINE_BI_CONSUMER);
    }
}