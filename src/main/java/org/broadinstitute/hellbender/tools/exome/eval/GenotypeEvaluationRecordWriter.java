package org.broadinstitute.hellbender.tools.exome.eval;

import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;

import java.io.File;
import java.io.IOException;
import java.util.stream.Collectors;

/**
 * Table file writer for {@link GenotypeEvaluationRecord} instances.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
class GenotypeEvaluationRecordWriter extends TableWriter<GenotypeEvaluationRecord> {

    private static final String SAMPLE_NAME_COLUMN = "SAMPLE";
    private static final String EVALUATION_CLASS_COLUMN = "CLASS";
    private static final String CONTIG_COLUMN = "CONTIG";
    private static final String START_COLUMN = "START";
    private static final String END_COLUMN = "END";
    private static final String TARGET_COUNT_COLUMN = "NTARGETS";
    private static final String FILTER_COLUMN = "FILTER";
    private static final String TRUTH_COLUMN = "TRUTH";
    private static final String CALL_COLUMN = "CALL";
    private static final String OUTCOME_COLUMN = "RESULT";

    private static final TableColumnCollection COLUMNS = new TableColumnCollection(
            CONTIG_COLUMN,
            START_COLUMN,
            END_COLUMN,
            TARGET_COUNT_COLUMN,
            SAMPLE_NAME_COLUMN,
            EVALUATION_CLASS_COLUMN,
            FILTER_COLUMN,
            OUTCOME_COLUMN,
            TRUTH_COLUMN,
            CALL_COLUMN
    );

    private static final String SEGMENT_SEPARATOR = ";";
    private static final String ATTRIBUTE_SEPARATOR = ":";
    private static final String NO_VALUE = ".";

    public GenotypeEvaluationRecordWriter(final File file) throws IOException {
        super(file, COLUMNS);
        writeComment("Possible classes in the CLASS column: ");
        for (final EvaluationClass clazz : EvaluationClass.values()) {
            writeComment(String.format("    %s : %s", clazz.acronym, clazz.name().replace("_", " ")));
        }
        writeComment(String.format("Content of %s and %s columns:", TRUTH_COLUMN, CALL_COLUMN));
        writeComment("   They may contain multiple segment information or none (represented with '.')");
        writeComment(String.format("   Each segment record is separated with a '%s' and each field in each segment is " +
                                      "separated by '%s'", SEGMENT_SEPARATOR, ATTRIBUTE_SEPARATOR));
        writeComment("   The segments attributes are in this order:");
        writeComment("       " + String.join(ATTRIBUTE_SEPARATOR, "Start", "End", "Call", "TargetCount", "Mean Cov.",
                "Std Cov.", "Some Quality", "Start Quality", "End Quality"));
    }

    @Override
    protected void composeLine(final GenotypeEvaluationRecord record,
                               final DataLine dataLine) {
        dataLine.append(record.interval.getContig())
                .append(record.interval.getStart())
                .append(record.interval.getEnd())
                .append(record.targetCount)
                .append(record.sample)
                .append(record.evaluationClass == null ? NO_VALUE : record.evaluationClass.toString())
                .append(record.getFilterString())
                .append(record.result.toString())
                .append(record.truth == null ? NO_VALUE : record.truth.toString())
                .append(record.calls.size() == 0 ? NO_VALUE : record.calls.stream()
                        .map(Object::toString)
                        .collect(Collectors.joining("; ")));
    }
}
