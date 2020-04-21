package org.broadinstitute.hellbender.tools.walkers.contamination;

import java.nio.file.Path;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
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
    private String sample;
    private double contamination;
    private double error;

    public ContaminationRecord(final String sample, final double contamination, final double error) {
        this.sample = sample;
        this.contamination = contamination;
        this.error = error;
    }

    public String getSample() { return sample; }

    public double getContamination() {
        return contamination;
    }

    public double getError() {
        return error;
    }

    //----- The following two public static methods read and write contamination files
    public static void writeToFile(final List<ContaminationRecord> records, final File outputTable) {
        try ( ContaminationRecord.ContaminationTableWriter writer = new ContaminationRecord.ContaminationTableWriter(IOUtils.toPath(outputTable)) ) {
            writer.writeAllRecords(records);
        } catch (IOException e){
            throw new UserException(String.format("Encountered an IO exception while writing to %s.", outputTable));
        }
    }

    public static List<ContaminationRecord> readFromFile(final File tableFile) {
        try( ContaminationTableReader reader = new ContaminationTableReader(IOUtils.toPath(tableFile)) ) {
            return reader.toList();
        } catch (IOException e){
            throw new UserException(String.format("Encountered an IO exception while reading from %s.", tableFile));
        }
    }

    //-------- The following methods are boilerplate for reading and writing contamination tables
    private static class ContaminationTableWriter extends TableWriter<ContaminationRecord> {
        private ContaminationTableWriter(final Path output) throws IOException {
            super(output, ContaminationTableColumn.COLUMNS);
        }

        @Override
        protected void composeLine(final ContaminationRecord record, final DataLine dataLine) {
            dataLine.set(ContaminationTableColumn.SAMPLE.toString(), record.getSample())
                .set(ContaminationTableColumn.CONTAMINATION.toString(), record.getContamination())
                    .set(ContaminationTableColumn.ERROR.toString(), record.getError());
        }
    }

    private static class ContaminationTableReader extends TableReader<ContaminationRecord> {
        public ContaminationTableReader(final Path path) throws IOException {
            super(path);
        }

        @Override
        protected ContaminationRecord createRecord(final DataLine dataLine) {
            final String sample = dataLine.get(ContaminationTableColumn.SAMPLE);
            final double contamination = dataLine.getDouble(ContaminationTableColumn.CONTAMINATION);
            final double error = dataLine.getDouble(ContaminationTableColumn.ERROR);
            return new ContaminationRecord(sample, contamination, error);
        }
    }

    private enum ContaminationTableColumn {
        SAMPLE("sample"),
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

        public static final TableColumnCollection COLUMNS = new TableColumnCollection(SAMPLE, CONTAMINATION, ERROR);
    }
}
