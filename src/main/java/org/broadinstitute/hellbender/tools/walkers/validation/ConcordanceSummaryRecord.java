package org.broadinstitute.hellbender.tools.walkers.validation;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableReader;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;

import java.io.File;
import java.io.IOException;

/**
 * Created by tsato on 2/8/17.
 */

public class ConcordanceSummaryRecord {
    private static final String VARIANT_TYPE_COLUMN_NAME = "type";
    private static final String TRUE_POSITIVE_COLUMN_NAME = "true-positive";
    private static final String FALSE_POSITIVE_COLUMN_NAME = "false-positive";
    private static final String FALSE_NEGATIVE_COLUMN_NAME = "false-negative";
    private static final String SENSITIVITY_COLUMN_NAME = "sensitivity";
    private static final String PRECISION_COLUMN_NAME = "precision";
    private static final String[] SUMMARY_TABLE_COLUMN_HEADER =
            {VARIANT_TYPE_COLUMN_NAME, TRUE_POSITIVE_COLUMN_NAME, FALSE_POSITIVE_COLUMN_NAME,
                    FALSE_NEGATIVE_COLUMN_NAME, SENSITIVITY_COLUMN_NAME, PRECISION_COLUMN_NAME};

    final VariantContext.Type type;
    final long truePositives;
    final long falsePositives;
    final long falseNegatives;

    public ConcordanceSummaryRecord(final VariantContext.Type type, final long truePositives, final long falsePositives, final long falseNegatives){
        this.type = type;
        this.truePositives = truePositives;
        this.falsePositives = falsePositives;
        this.falseNegatives = falseNegatives;
    }

    public VariantContext.Type getVariantType() { return type; }

    public long getTruePositives() { return truePositives; }

    public long getFalsePositives() { return falsePositives; }

    public long getFalseNegatives() { return falseNegatives; }

    public double getSensitivity() { return (double) truePositives / (truePositives + falseNegatives); }

    public double getPrecision() { return (double) truePositives / (truePositives + falsePositives); }

    public static class Writer extends TableWriter<ConcordanceSummaryRecord> {
        private Writer(final File output) throws IOException {
            super(output, new TableColumnCollection(SUMMARY_TABLE_COLUMN_HEADER));
        }

        @Override
        protected void composeLine(final ConcordanceSummaryRecord record, final DataLine dataLine) {
            dataLine.set(VARIANT_TYPE_COLUMN_NAME, record.getVariantType().toString())
                    .set(TRUE_POSITIVE_COLUMN_NAME, record.getTruePositives())
                    .set(FALSE_POSITIVE_COLUMN_NAME, record.getFalsePositives())
                    .set(FALSE_NEGATIVE_COLUMN_NAME, record.getFalseNegatives())
                    .set(SENSITIVITY_COLUMN_NAME, record.getSensitivity())
                    .set(PRECISION_COLUMN_NAME, record.getPrecision());
        }
    }

    public static Writer getWriter(final File outputTable){
        try {
            Writer writer = new Writer(outputTable);
            return writer;
        } catch (IOException e){
            throw new UserException(String.format("Encountered an IO exception while reading from %s.", outputTable), e);
        }
    }

    public static class Reader extends TableReader<ConcordanceSummaryRecord> {
        public Reader(final File summary) throws IOException {
            super(summary);
        }

        @Override
        protected ConcordanceSummaryRecord createRecord(final DataLine dataLine) {
            final VariantContext.Type type = VariantContext.Type.valueOf(dataLine.get(VARIANT_TYPE_COLUMN_NAME));
            final long truePositives = Long.parseLong(dataLine.get(TRUE_POSITIVE_COLUMN_NAME));
            final long falsePositives = Long.parseLong(dataLine.get(FALSE_POSITIVE_COLUMN_NAME));
            final long falseNegatives = Long.parseLong(dataLine.get(FALSE_NEGATIVE_COLUMN_NAME));

            return new ConcordanceSummaryRecord(type, truePositives, falsePositives, falseNegatives);
        }
    }
}

