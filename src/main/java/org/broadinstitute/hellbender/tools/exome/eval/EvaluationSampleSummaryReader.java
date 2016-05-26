package org.broadinstitute.hellbender.tools.exome.eval;

import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableReader;

import java.io.File;
import java.io.IOException;

/**
 * Table reader for {@link EvaluationSampleSummaryRecord} instances.
 */
public final class EvaluationSampleSummaryReader extends TableReader<EvaluationSampleSummaryRecord> {

    public EvaluationSampleSummaryReader(final File file) throws IOException {
        super(file);
    }

    @Override
    protected EvaluationSampleSummaryRecord createRecord(final DataLine dataLine) {
        final String sample = dataLine.get(EvaluationSampleSummaryWriter.SAMPLE_NAME_COLUMN_NAME);
        final EvaluationSampleSummaryRecord result = new EvaluationSampleSummaryRecord(sample);
        for (final EvaluationClass evaluationClass : EvaluationClass.values()) {
            result.set(evaluationClass, dataLine.getInt(evaluationClass.acronym));
        }
        return result;
    }
}
