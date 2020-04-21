package org.broadinstitute.hellbender.tools.walkers.contamination;

import htsjdk.samtools.util.Locatable;
import java.nio.file.Path;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.tsv.*;

import java.io.File;
import java.io.IOException;
import java.util.List;

/**
 * Created by David Benjamin on 2/13/17.
 */
public class MinorAlleleFractionRecord implements Locatable {
    private final SimpleInterval segment;
    private double minorAlleleFraction;

    public MinorAlleleFractionRecord(final SimpleInterval segment, final double minorAlleleFraction) {
        this.segment = segment;
        this.minorAlleleFraction = minorAlleleFraction;
    }

    @Override
    public String getContig() { return segment.getContig(); }

    @Override
    public int getStart() { return segment.getStart(); }

    @Override
    public int getEnd() { return segment.getEnd(); }

    public SimpleInterval getSegment() {
        return segment;
    }

    public double getMinorAlleleFraction() {
        return minorAlleleFraction;
    }

    //----- The following two public static methods read and write contamination files
    public static void writeToFile(final String sample, final List<MinorAlleleFractionRecord> records, final File outputTable) {
        try ( MinorAlleleFractionTableWriter writer = new MinorAlleleFractionTableWriter(IOUtils.toPath(outputTable)) ) {
            writer.writeMetadata(TableUtils.SAMPLE_METADATA_TAG, sample);
            writer.writeAllRecords(records);
        } catch (IOException e){
            throw new UserException(String.format("Encountered an IO exception while writing to %s.", outputTable));
        }
    }

    public static ImmutablePair<String, List<MinorAlleleFractionRecord>> readFromFile(final File tableFile) {
        return readFromPath(IOUtils.toPath(tableFile));
    }

    public static ImmutablePair<String, List<MinorAlleleFractionRecord>> readFromPath(final Path tablePath) {
        try( MinorAlleleFractionTableReader reader = new MinorAlleleFractionTableReader(tablePath) ) {
            final List<MinorAlleleFractionRecord> list = reader.toList();
            return ImmutablePair.of(reader.getMetadata().get(TableUtils.SAMPLE_METADATA_TAG), list);
        } catch (IOException e){
            throw new UserException(String.format("Encountered an IO exception while reading from %s.", tablePath));
        }
    }

    //-------- The following methods are boilerplate for reading and writing contamination tables
    private static class MinorAlleleFractionTableWriter extends TableWriter<MinorAlleleFractionRecord> {
        private MinorAlleleFractionTableWriter(final Path output) throws IOException {
            super(output, MinorAlleleFractionTableColumn.COLUMNS);
        }

        @Override
        protected void composeLine(final MinorAlleleFractionRecord record, final DataLine dataLine) {
            final SimpleInterval segment = record.getSegment();
            dataLine.set(MinorAlleleFractionTableColumn.CONTIG.toString(), segment.getContig())
                    .set(MinorAlleleFractionTableColumn.START.toString(), segment.getStart())
                    .set(MinorAlleleFractionTableColumn.END.toString(), segment.getEnd())
                    .set(MinorAlleleFractionTableColumn.MAF.toString(), record.getMinorAlleleFraction());
        }
    }

    public static class MinorAlleleFractionTableReader extends TableReader<MinorAlleleFractionRecord> {
        public MinorAlleleFractionTableReader(final Path path) throws IOException {
            super(path);
        }

        @Override
        protected MinorAlleleFractionRecord createRecord(final DataLine dataLine) {
            final String contig = dataLine.get(MinorAlleleFractionTableColumn.CONTIG);
            final int start = dataLine.getInt(MinorAlleleFractionTableColumn.START);
            final int end = dataLine.getInt(MinorAlleleFractionTableColumn.END);
            final double maf = dataLine.getDouble(MinorAlleleFractionTableColumn.MAF);
            return new MinorAlleleFractionRecord(new SimpleInterval(contig, start, end), maf);
        }
    }

    private enum MinorAlleleFractionTableColumn {
        CONTIG("contig"),
        START("start"),
        END("end"),
        MAF("minor_allele_fraction");

        private final String columnName;

        MinorAlleleFractionTableColumn(final String columnName) {
            this.columnName = Utils.nonNull(columnName);
        }

        @Override
        public String toString() {
            return columnName;
        }

        public static final TableColumnCollection COLUMNS = new TableColumnCollection(CONTIG, START, END, MAF);
    }
}
