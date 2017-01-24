package org.broadinstitute.hellbender.tools.copynumber;

import org.broadinstitute.hellbender.tools.copynumber.coverage.readcount.ReadCountData;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;

import java.io.File;
import java.io.IOException;
import java.io.Writer;

/**
 * Created by asmirnov on 7/24/17.
 */
public class BinnedReadCountWriter extends TableWriter<ReadCountData> {

    public BinnedReadCountWriter(final File file, final TableColumnCollection columns) throws IOException {
        super(file, columns);
    }

    public BinnedReadCountWriter(final Writer writer, final TableColumnCollection columns) throws IOException {
        super(writer, columns);
    }

    @Override
    protected void composeLine(ReadCountData record, DataLine dataLine) {
        SimpleInterval interval = record.getInterval();
        dataLine.append(interval.getContig())
                .append(interval.getStart())
                .append(interval.getEnd());
        record.appendCountsTo(dataLine);
    }
}
