package org.broadinstitute.hellbender.utils.hmm.segmentation;

import org.broadinstitute.hdf5.Utils;
import org.broadinstitute.hellbender.tools.exome.SegmentTableColumn;
import org.broadinstitute.hellbender.tools.exome.Target;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableReader;

import java.io.File;
import java.io.IOException;
import java.util.function.Function;

/**
 * Reads {@link HiddenStateSegmentRecord} instances from a tab separated table file.
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public class HiddenStateSegmentRecordReader<S, T extends Target> extends TableReader<HiddenStateSegmentRecord<S, T>> {

    /**
     * The inverse call string factory maps a call string to its corresponding hidden state
     */
    private final Function<String, S> callStringParser;

    /**
     * Opens a reader on an a preexisting file.
     *
     * @param file the source file where to read from.
     * @param callStringParser  a map from call string to a hidden state
     *
     * @throws IllegalArgumentException if {@code file} is {@code null}.
     * @throws IOException if there is an issue trying to read the contents of the file.
     * @throws RuntimeException if there is a formatting issue within the file.
     */
    public HiddenStateSegmentRecordReader(final File file,
                                          final Function<String, S> callStringParser) throws IOException {
        super(file);
        this.callStringParser = Utils.nonNull(callStringParser);
    }

    @Override
    protected HiddenStateSegmentRecord<S, T> createRecord(DataLine dataLine) {
        final SimpleInterval interval = new SimpleInterval(
                dataLine.get(SegmentTableColumn.CONTIG),
                dataLine.getInt(SegmentTableColumn.START),
                dataLine.getInt(SegmentTableColumn.END)
        );
        final HiddenStateSegment<S, T> segment = new HiddenStateSegment<>(
                interval, dataLine.getInt(SegmentTableColumn.NUM_TARGETS),
                dataLine.getDouble(SegmentTableColumn.MEAN),
                dataLine.getDouble(SegmentTableColumn.SD),
                callStringParser.apply(dataLine.get(SegmentTableColumn.CALL)),
                dataLine.getDouble(SegmentTableColumn.EXACT_QUALITY),
                dataLine.getDouble(SegmentTableColumn.SOME_QUALITY),
                dataLine.getDouble(SegmentTableColumn.START_QUALITY),
                dataLine.getDouble(SegmentTableColumn.END_QUALITY),
                dataLine.getDouble(SegmentTableColumn.EVENT_QUALITY));

        return new HiddenStateSegmentRecord<>(
                dataLine.get(SegmentTableColumn.SAMPLE), segment);
    }
}
