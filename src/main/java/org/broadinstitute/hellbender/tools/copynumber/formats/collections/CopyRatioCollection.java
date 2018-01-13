package org.broadinstitute.hellbender.tools.copynumber.formats.collections;

import htsjdk.samtools.util.OverlapDetector;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SampleLocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.CopyRatio;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;

import java.io.File;
import java.util.List;
import java.util.function.BiConsumer;
import java.util.function.Function;
import java.util.stream.Collectors;

/**
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class CopyRatioCollection extends AbstractSampleLocatableCollection<CopyRatio> {
    //note to developers: repeat the column headers in Javadoc so that they are viewable when linked
    /**
     * CONTIG, START, END, LOG2_COPY_RATIO
     */
    enum CopyRatioTableColumn {
        CONTIG,
        START,
        END,
        LOG2_COPY_RATIO;

        static final TableColumnCollection COLUMNS = new TableColumnCollection((Object[]) values());
    }

    private static final Function<DataLine, CopyRatio> COPY_RATIO_RECORD_FROM_DATA_LINE_DECODER = dataLine -> {
        final String contig = dataLine.get(CopyRatioTableColumn.CONTIG);
        final int start = dataLine.getInt(CopyRatioTableColumn.START);
        final int end = dataLine.getInt(CopyRatioTableColumn.END);
        final double copyRatio = dataLine.getDouble(CopyRatioTableColumn.LOG2_COPY_RATIO);
        final SimpleInterval interval = new SimpleInterval(contig, start, end);
        return new CopyRatio(interval, copyRatio);
    };

    private static final BiConsumer<CopyRatio, DataLine> COPY_RATIO_RECORD_TO_DATA_LINE_ENCODER = (copyRatio, dataLine) ->
        dataLine.append(copyRatio.getInterval().getContig())
                .append(copyRatio.getInterval().getStart())
                .append(copyRatio.getInterval().getEnd())
                .append(formatDouble(copyRatio.getLog2CopyRatioValue()));

    public CopyRatioCollection(final File inputFile) {
        super(inputFile, CopyRatioTableColumn.COLUMNS, COPY_RATIO_RECORD_FROM_DATA_LINE_DECODER, COPY_RATIO_RECORD_TO_DATA_LINE_ENCODER);
    }

    public CopyRatioCollection(final SampleLocatableMetadata metadata,
                               final List<CopyRatio> copyRatios) {
        super(metadata, copyRatios, CopyRatioTableColumn.COLUMNS, COPY_RATIO_RECORD_FROM_DATA_LINE_DECODER, COPY_RATIO_RECORD_TO_DATA_LINE_ENCODER);
    }

    public List<Double> getLog2CopyRatioValues() {
        return getRecords().stream().map(CopyRatio::getLog2CopyRatioValue).collect(Collectors.toList());
    }

    /**
     * The midpoint is used to characterize the interval for the purposes of determining overlaps
     * so that each copy-ratio interval will be uniquely contained in a single segment.
     */
    public OverlapDetector<CopyRatio> getMidpointOverlapDetector() {
        return OverlapDetector.create(getRecords().stream()
                .map(cr -> new CopyRatio(cr.getMidpoint(), cr.getLog2CopyRatioValue()))
                .collect(Collectors.toList()));
    }
}
