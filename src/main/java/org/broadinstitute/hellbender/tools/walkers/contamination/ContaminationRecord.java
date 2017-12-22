package org.broadinstitute.hellbender.tools.walkers.contamination;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableReader;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;

import java.io.File;
import java.io.IOException;
import java.util.List;

/**
 * Created by David Benjamin on 2/13/17.
 */
public class ContaminationRecord {
    private String level;
    private double contamination;
    private double error;

    public ContaminationRecord(final String level, final double contamination, final double error) {
        this.level = level;
        this.contamination = contamination;
        this.error = error;
    }

    public String getLevel() { return level; }

    public double getContamination() {
        return contamination;
    }

    public double getError() {
        return error;
    }

    //----- The following two public static methods read and write contamination files
    public static void writeToFile(final List<ContaminationRecord> records, final File outputTable) {
        try ( ContaminationRecord.ContaminationTableWriter writer = new ContaminationRecord.ContaminationTableWriter(outputTable) ) {
            writer.writeAllRecords(records);
        } catch (IOException e){
            throw new UserException(String.format("Encountered an IO exception while writing to %s.", outputTable));
        }
    }

    public static List<ContaminationRecord> readFromFile(final File tableFile) {
        try( ContaminationTableReader reader = new ContaminationTableReader(tableFile) ) {
            return reader.toList();
        } catch (IOException e){
            throw new UserException(String.format("Encountered an IO exception while reading from %s.", tableFile));
        }
    }

    //-------- The following methods are boilerplate for reading and writing contamination tables
    private static class ContaminationTableWriter extends TableWriter<ContaminationRecord> {
        private ContaminationTableWriter(final File output) throws IOException {
            super(output, ContaminationTableColumn.COLUMNS);
        }

        @Override
        protected void composeLine(final ContaminationRecord record, final DataLine dataLine) {
            dataLine.set(ContaminationTableColumn.LEVEL.toString(), record.getLevel())
                    .set(ContaminationTableColumn.CONTAMINATION.toString(), record.getContamination())
                    .set(ContaminationTableColumn.ERROR.toString(), record.getError());
        }
    }

    private static class ContaminationTableReader extends TableReader<ContaminationRecord> {
        public ContaminationTableReader(final File file) throws IOException {
            super(file);
        }

        @Override
        protected ContaminationRecord createRecord(final DataLine dataLine) {
            final String sample = dataLine.get(ContaminationTableColumn.LEVEL);
            final double contamination = dataLine.getDouble(ContaminationTableColumn.CONTAMINATION);
            final double error = dataLine.getDouble(ContaminationTableColumn.ERROR);
            return new ContaminationRecord(sample, contamination, error);
        }
    }

    private enum ContaminationTableColumn {
        LEVEL("level"),
        CONTAMINATION("contamination"),
        ERROR("error");

        private final String columnName;

        ContaminationTableColumn(final String columnName) {
            this.columnName = Utils.nonNull(columnName);
        }

        @Override
        public String toString() {
            return columnName;
        }

        public static final TableColumnCollection COLUMNS = new TableColumnCollection(LEVEL, CONTAMINATION, ERROR);
    }

    public enum Level {
        WHOLE_BAM("whole_bam"),
        READ_GROUP("read_group"),
        LANE("lane");

        private String value;

        Level(final String value) { this.value = value; }

        @Override
        public String toString() { return value; }
    }
}
