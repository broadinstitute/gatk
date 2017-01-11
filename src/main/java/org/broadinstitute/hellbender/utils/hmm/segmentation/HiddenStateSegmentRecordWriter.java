package org.broadinstitute.hellbender.utils.hmm.segmentation;

import org.broadinstitute.hellbender.tools.exome.SegmentTableColumn;
import org.broadinstitute.hellbender.tools.exome.Target;
import org.broadinstitute.hellbender.utils.hmm.interfaces.CallStringProducer;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;

import java.io.File;
import java.io.IOException;

/**
 * Table writer for {@link HiddenStateSegmentRecord}, thus {@link HiddenStateSegment}
 * by extension.
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public class HiddenStateSegmentRecordWriter<S extends CallStringProducer, T extends Target> extends TableWriter<HiddenStateSegmentRecord<S, T>> {

    public HiddenStateSegmentRecordWriter(final File file) throws IOException {
        super(file, SegmentTableColumn.GERMLINE_CALL_COLUMNS);
    }

    @Override
    protected void composeLine(final HiddenStateSegmentRecord<S, T> record, final DataLine dataLine) {
        final HiddenStateSegment<S, T> segment = record.getSegment();
        dataLine.append(record.getSampleName())
                .append(segment.getContig())
                .append(segment.getStart())
                .append(segment.getEnd())
                .append(segment.getCall().getCallString())
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
