package org.broadinstitute.hellbender.tools.exome.titanconversion;

import org.broadinstitute.hellbender.tools.exome.ReadCountRecord;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;

import java.io.File;
import java.io.IOException;

public class TitanCopyRatioEstimateWriter extends TableWriter<ReadCountRecord> {

    public TitanCopyRatioEstimateWriter(final File file) throws IOException {
        super(file, new TableColumnCollection(TitanCopyRatioEstimateColumns.FULL_COLUMN_NAME_ARRAY));
    }

    @Override
    protected void composeLine(ReadCountRecord record, DataLine dataLine) {

        // chr	start	end	log2_TNratio_corrected
        dataLine.append(record.getContig())
                .append(record.getStart())
                .append(record.getEnd())
                .append(record.getDoubleCounts()[0]);
    }
}
