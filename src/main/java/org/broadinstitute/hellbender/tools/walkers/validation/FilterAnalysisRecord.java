package org.broadinstitute.hellbender.tools.walkers.validation;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableReader;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.List;

public class FilterAnalysisRecord {
    private final String filter;
    private int trueNegativeCount;
    private int falseNegativeCount;
    private int uniqueTrueNegativeCount;
    private int uniqueFalseNegativeCount;

    public FilterAnalysisRecord(final String filter, final int trueNegativeCount, final int falseNegativeCount, final int uniqueTrueNegativeCount, final int uniqueFalseNegativeCount) {
        this.filter = filter;
        this.trueNegativeCount = trueNegativeCount;
        this.falseNegativeCount = falseNegativeCount;
        this.uniqueTrueNegativeCount = uniqueTrueNegativeCount;
        this.uniqueFalseNegativeCount = uniqueFalseNegativeCount;
    }

    public String getFilter() { return filter; }

    public int getTrueNegativeCount() { return trueNegativeCount; }

    public int getFalseNegativeCount() { return falseNegativeCount; }

    public int getUniqueTrueNegativeCount() { return uniqueTrueNegativeCount; }

    public int getUniqueFalseNegativeCount() { return uniqueFalseNegativeCount; }

    public void incrementTrueNegative() { trueNegativeCount++; }
    public void incrementFalseNegative() { falseNegativeCount++; }
    public void incrementUniqueTrueNegative() { uniqueTrueNegativeCount++; }
    public void incrementUniqueFalseNegative() { uniqueFalseNegativeCount++; }

    //----- The following two public static methods read and write contamination files
    public static void writeToFile(final Collection<FilterAnalysisRecord> records, final File outputTable) {
        try ( FilterAnalysisTableWriter writer = new FilterAnalysisTableWriter(outputTable) ) {
            writer.writeAllRecords(records);
        } catch (IOException e){
            throw new UserException(String.format("Encountered an IO exception while writing to %s.", outputTable));
        }
    }

    public static List<FilterAnalysisRecord> readFromFile(final File tableFile) {
        try( FilterAnalysisTableReader reader = new FilterAnalysisTableReader(tableFile) ) {
            return reader.toList();
        } catch (IOException e){
            throw new UserException(String.format("Encountered an IO exception while reading from %s.", tableFile));
        }
    }

    //-------- The following methods are boilerplate for reading and writing contamination tables
    private static class FilterAnalysisTableWriter extends TableWriter<FilterAnalysisRecord> {
        private FilterAnalysisTableWriter(final File output) throws IOException {
            super(output, FilterAnalysisTableColumn.COLUMNS);
        }

        @Override
        protected void composeLine(final FilterAnalysisRecord record, final DataLine dataLine) {
            dataLine.set(FilterAnalysisTableColumn.FILTER.toString(), record.getFilter())
                    .set(FilterAnalysisTableColumn.TRUE_NEGATIVES.toString(), record.getTrueNegativeCount())
                    .set(FilterAnalysisTableColumn.FALSE_NEGATIVES.toString(), record.getFalseNegativeCount())
                    .set(FilterAnalysisTableColumn.UNIQUE_TRUE_NEGATIVES.toString(), record.getUniqueTrueNegativeCount())
                    .set(FilterAnalysisTableColumn.UNIQUE_FALSE_NEGATIVES.toString(), record.getUniqueFalseNegativeCount());
        }
    }

    private static class FilterAnalysisTableReader extends TableReader<FilterAnalysisRecord> {
        public FilterAnalysisTableReader(final File file) throws IOException {
            super(file);
        }

        @Override
        protected FilterAnalysisRecord createRecord(final DataLine dataLine) {
            final String filter = dataLine.get(FilterAnalysisTableColumn.FILTER);
            final int trueNegatives = dataLine.getInt(FilterAnalysisTableColumn.TRUE_NEGATIVES);
            final int falseNegatives = dataLine.getInt(FilterAnalysisTableColumn.FALSE_NEGATIVES);
            final int uniqueTrueNegatives = dataLine.getInt(FilterAnalysisTableColumn.UNIQUE_TRUE_NEGATIVES);
            final int uniqueFalseNegatives = dataLine.getInt(FilterAnalysisTableColumn.UNIQUE_FALSE_NEGATIVES);
            return new FilterAnalysisRecord(filter, trueNegatives, falseNegatives, uniqueTrueNegatives, uniqueFalseNegatives);
        }
    }

    private enum FilterAnalysisTableColumn {
        FILTER("filter"),
        TRUE_NEGATIVES("tn"),
        FALSE_NEGATIVES("fn"),
        UNIQUE_TRUE_NEGATIVES("uniq_tn"),
        UNIQUE_FALSE_NEGATIVES("uniq_fn");

        private final String columnName;

        FilterAnalysisTableColumn(final String columnName) {
            this.columnName = Utils.nonNull(columnName);
        }

        @Override
        public String toString() {
            return columnName;
        }

        public static final TableColumnCollection COLUMNS = new TableColumnCollection(FILTER, TRUE_NEGATIVES, FALSE_NEGATIVES, UNIQUE_TRUE_NEGATIVES, UNIQUE_FALSE_NEGATIVES);
    }
}
