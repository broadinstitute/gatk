package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;

import java.io.File;
import java.io.IOException;

/**
 * Table writer for {@link CopyNumberTriStateSegmentRecord}, thus {@link CopyNumberTriStateSegment}
 * by extension.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class CopyNumberTriStateSegmentRecordWriter extends TableWriter<CopyNumberTriStateSegmentRecord> {

    public CopyNumberTriStateSegmentRecordWriter(final File file) throws IOException {
        super(file, SegmentTableColumn.GERMLINE_CALL_COLUMNS);
    }

    @Override
    protected void composeLine(final CopyNumberTriStateSegmentRecord record, final DataLine dataLine) {
        final CopyNumberTriStateSegment segment = record.getSegment();
        dataLine.append(record.getSampleName())
                .append(segment.getContig())
                .append(segment.getStart())
                .append(segment.getEnd())
                .append(segment.getCall().callString)
                .append(segment.getTargetCount())
                .append(formatDouble(segment.getMean()))
                .append(formatDouble(segment.getStdev()))
                .append(formatDouble(segment.getExactQuality()))
                .append(formatDouble(segment.getSomeQuality()))
                .append(formatDouble(segment.getStartQuality()))
                .append(formatDouble(segment.getEndQuality()))
                .append(formatDouble(segment.getEventQuality()));
    }

    private static String formatDouble(final double d) {
        return String.format("%.4f", d);
    }
}
