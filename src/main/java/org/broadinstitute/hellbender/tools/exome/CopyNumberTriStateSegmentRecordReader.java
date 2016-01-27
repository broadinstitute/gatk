package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.hmm.CopyNumberTriState;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableReader;

import java.io.File;
import java.io.IOException;

/**
 * Reads {@link CopyNumberTriStateSegmentRecord} instances from a tab separated table file.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public class CopyNumberTriStateSegmentRecordReader extends TableReader<CopyNumberTriStateSegmentRecord> {

    /**
     * Opens a reader on an a preexisting file.
     *
     * @param file the source file where to read from.
     *
     * @throws IllegalArgumentException if {@code file} is {@code null}.
     * @throws IOException if there is an issue trying to read the contents of the file.
     * @throws RuntimeException if there is a formatting issue within the file.
     */
    public CopyNumberTriStateSegmentRecordReader(final File file) throws IOException {
        super(file);
    }

    @Override
    protected CopyNumberTriStateSegmentRecord createRecord(DataLine dataLine) {
        final SimpleInterval interval = new SimpleInterval(
                dataLine.get(SegmentTableColumns.CONTIG),
                dataLine.getInt(SegmentTableColumns.START),
                dataLine.getInt(SegmentTableColumns.END)
        );
        final CopyNumberTriStateSegment segment = new CopyNumberTriStateSegment(
                interval, dataLine.getInt(SegmentTableColumns.NUM_TARGETS),
                dataLine.getDouble(SegmentTableColumns.MEAN),
                dataLine.getDouble(SegmentTableColumns.SD),
                CopyNumberTriState.fromCallString(dataLine.get(SegmentTableColumns.CALL)),
                dataLine.getDouble(SegmentTableColumns.EXACT_QUALITY),
                dataLine.getDouble(SegmentTableColumns.SOME_QUALITY),
                dataLine.getDouble(SegmentTableColumns.START_QUALITY),
                dataLine.getDouble(SegmentTableColumns.END_QUALITY),
                dataLine.getDouble(SegmentTableColumns.EVENT_QUALITY));

        return new CopyNumberTriStateSegmentRecord(
                dataLine.get(SegmentTableColumns.SAMPLE), segment);
    }
}
