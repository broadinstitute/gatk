package org.broadinstitute.hellbender.tools.walkers.validation;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableReader;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.List;

/**
 * Simple class for storing a sample and its mixing fraction within a pooled bam.
 *
 * Created by David Benjamin on 1/31/17.
 */
public class MixingFraction {
    private final String sample;
    private final double mixingFraction;

    public MixingFraction(final String sample, final double mixingFraction) {
        this.sample = sample;
        this.mixingFraction = mixingFraction;
    }

    public String getSample() {
        return sample;
    }

    public double getMixingFraction() {
        return mixingFraction;
    }

    public static List<MixingFraction> readMixingFractions(final File file) {
        Utils.regularReadableUserFile(file);
        try (final MixingFractionReader reader = new MixingFractionReader(file)) {
            return reader.toList();
        } catch (final FileNotFoundException ex) {
            throw new UserException.CouldNotReadInputFile("Mixing fraction table file not found.", ex);
        } catch (final IOException ex) {
            throw new UserException.CouldNotReadInputFile("Could not read input table.", ex);
        }
    }

    public static void writeMixingFractions(List<MixingFraction> mixingFractions, final File file) {
        try (MixingFractionWriter writer = new MixingFractionWriter(file)) {
            writer.writeAllRecords(mixingFractions);
        } catch (IOException e){
            throw new UserException(String.format("Encountered an IO exception while trying to create output file %s.", file));
        }
    }

    private enum MixingFractionTableColumn {
        SAMPLE("SAMPLE"),
        MIXING_FRACTION("MIXING_FRACTION");

        private final String columnName;  //store the column names

        MixingFractionTableColumn(final String columnName) {
            this.columnName = Utils.nonNull(columnName);
        }

        @Override
        public String toString() {
            return columnName;
        }

        public static final TableColumnCollection COLUMNS = new TableColumnCollection((Object[]) values());
    }

    private static class MixingFractionWriter extends TableWriter<MixingFraction> {
        public MixingFractionWriter(final File output) throws IOException {
            super(output, MixingFractionTableColumn.COLUMNS);
        }

        @Override
        protected void composeLine(final MixingFraction record, final DataLine dataLine) {
            dataLine.set(MixingFractionTableColumn.SAMPLE.ordinal(), record.getSample())
                    .set(MixingFractionTableColumn.MIXING_FRACTION.ordinal(), record.getMixingFraction());
        }
    }

    private static class MixingFractionReader extends TableReader<MixingFraction> {
        public MixingFractionReader(final File file) throws IOException {
            super(file);
        }

        @Override
        protected MixingFraction createRecord(final DataLine dataLine) {
            final String sample = dataLine.get(MixingFractionTableColumn.SAMPLE);
            final double mixingFraction = dataLine.getDouble(MixingFractionTableColumn.MIXING_FRACTION);
            return new MixingFraction(sample, mixingFraction);
        }
    }
}
