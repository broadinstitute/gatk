package org.broadinstitute.hellbender.tools.walkers.mutect;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableReader;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;

import java.io.File;
import java.io.IOException;
import java.util.List;

public class MutectStats {
    private String statistic;
    private double value;

    public MutectStats(final String statistic, final double value) {
        this.statistic = statistic;
        this.value = value;
    }

    public String getStatistic() { return statistic; }

    public double getValue() { return value; }
    

    //----- The following two public static methods read and write contamination files
    public static void writeToFile(final List<MutectStats> records, final File outputTable) {
        try ( MutectStats.MutectStatsWriter writer = new MutectStats.MutectStatsWriter(outputTable) ) {
            writer.writeAllRecords(records);
        } catch (IOException e){
            throw new UserException(String.format("Encountered an IO exception while writing to %s.", outputTable));
        }
    }

    public static List<MutectStats> readFromFile(final File tableFile) {
        try( MutectStats.MutectStatsReader reader = new MutectStats.MutectStatsReader(tableFile) ) {
            return reader.toList();
        } catch (IOException e){
            throw new UserException(String.format("Encountered an IO exception while reading from %s.", tableFile));
        }
    }

    //-------- The following methods are boilerplate for reading and writing contamination tables
    private static class MutectStatsWriter extends TableWriter<MutectStats> {
        private MutectStatsWriter(final File output) throws IOException {
            super(output, MutectStats.MutectStatsColumn.COLUMNS);
        }

        @Override
        protected void composeLine(final MutectStats record, final DataLine dataLine) {
            dataLine.set(MutectStats.MutectStatsColumn.STATISTIC.toString(), record.getStatistic())
                    .set(MutectStats.MutectStatsColumn.VALUE.toString(), record.getValue());
        }
    }

    private static class MutectStatsReader extends TableReader<MutectStats> {
        public MutectStatsReader(final File file) throws IOException {
            super(file);
        }

        @Override
        protected MutectStats createRecord(final DataLine dataLine) {
            final String sample = dataLine.get(MutectStats.MutectStatsColumn.STATISTIC);
            final double contamination = dataLine.getDouble(MutectStats.MutectStatsColumn.VALUE);
            return new MutectStats(sample, contamination);
        }
    }

    private enum MutectStatsColumn {
        STATISTIC("statistic"),
        VALUE("value");

        private final String columnName;

        MutectStatsColumn(final String columnName) {
            this.columnName = Utils.nonNull(columnName);
        }

        @Override
        public String toString() {
            return columnName;
        }

        public static final TableColumnCollection COLUMNS = new TableColumnCollection(STATISTIC, VALUE);
    }
}
