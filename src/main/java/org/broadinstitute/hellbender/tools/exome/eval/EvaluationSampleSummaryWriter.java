package org.broadinstitute.hellbender.tools.exome.eval;

import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * Evaluation class per sample record writer.
 */
public final class EvaluationSampleSummaryWriter extends TableWriter<EvaluationSampleSummaryRecord> {

    public static final String SAMPLE_NAME_COLUMN_NAME = "SAMPLE";
    public static final String TOTAL_COLUMN_NAME = "TOTAL";

    private final EvaluationSampleSummaryRecord overallRecord;

    private static final TableColumnCollection COLUMNS;

    static {
        final List<String> columnNames = new ArrayList<>(EvaluationClass.values().length + 2);
        columnNames.add(SAMPLE_NAME_COLUMN_NAME);
        for (final EvaluationClass evalClass : EvaluationClass.values()) {
            columnNames.add(evalClass.acronym);
        }
        columnNames.add(TOTAL_COLUMN_NAME);
        COLUMNS = new TableColumnCollection(columnNames);
    }

    EvaluationSampleSummaryWriter(final File file, final String overallSample) throws IOException {
        super(file, COLUMNS);
        this.overallRecord = new EvaluationSampleSummaryRecord(overallSample);
    }

    @Override
    public void writeRecord(final EvaluationSampleSummaryRecord record) throws IOException {
        overallRecord.aggregate(record);
        super.writeRecord(record);
    }

    @Override
    protected void composeLine(final EvaluationSampleSummaryRecord record, final DataLine dataLine) {
        dataLine.append(record.getSample());
        for (final EvaluationClass evalClass : EvaluationClass.values()) {
            dataLine.append(record.get(evalClass));
        }
        dataLine.append(record.getTotal());
    }

    public void writeOverallRecord() throws IOException {
        super.writeRecord(overallRecord);
    }
}
