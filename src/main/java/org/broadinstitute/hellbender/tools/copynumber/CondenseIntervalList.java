package org.broadinstitute.hellbender.tools.copynumber;

import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.IOUtil;
import org.apache.commons.lang.math.NumberUtils;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableReader;
import org.broadinstitute.hellbender.utils.tsv.TableUtils;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;
import picard.cmdline.programgroups.IntervalsManipulationProgramGroup;

import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;

/**
 * Merges consecutive intervals in a interval_list formatted file
 * that happen to overlap or be adjacent.
 */
@CommandLineProgramProperties(summary = "merges adjacent or overlapping consecutive intervals in .interval_list formatted files",
        oneLineSummary = "merges adjacent or overlapping consecutive intervals in .interval_list formatted files" ,
        programGroup = IntervalsManipulationProgramGroup.class)
@BetaFeature
public class CondenseIntervalList extends CommandLineProgram {

    public static final String CONDENSATION_FACTOR_LONG_NAME = "condensation-factor";
    public static final String MIN_OUT_INTERVAL_LENGTH_LONG_NAME = "min-out-interval-length";

    private enum IntervalsListColumns {
        CONTIG,
        START,
        END,
        STRAND,
        TARGET
    }

    @Argument(
            doc = "Path to the input counts file",
            fullName = StandardArgumentDefinitions.INPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME
    )
    private String input;

    @Argument(
            doc = "Path to the output counts file",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME
    )
    private String output;

    @Argument(
            doc = "Indicates how many input read count record must be condense into each output record. By default it condenses as many as possible",
            fullName = CONDENSATION_FACTOR_LONG_NAME,
            shortName = CONDENSATION_FACTOR_LONG_NAME,
            minValue = 1,
            optional = true
    )
    private int condensationFactor = Integer.MAX_VALUE;

    @Argument(
            doc = "expected bin-length after condensation, those output bin that do not match this length will be excluded. " +
                    "The default 0 means that such constrain is not be applied",
            fullName = MIN_OUT_INTERVAL_LENGTH_LONG_NAME,
            minValue = 1,
            optional = true
    )
    private int minimumOutputBinLength = 1;

    @Override
    public String doWork() {

        try (final IntervalsListReader reader = new IntervalsListReader(input);
             final IntervalsListWriter writer = new IntervalsListWriter(output, reader.columns())) {
            IntervalsListRecord record = reader.readRecord();
            final CondensationPan pan = CondensationPan.start(condensationFactor, reader.columns(), record);
            writer.addHeaderLines(reader.headerLines());
            writer.mustOutputColumnNameLine(reader.columnNameLineIsPresent());
            while ((record = reader.readRecord()) != null) {
                if (!pan.add(record)) {
                    final IntervalsListRecord condensate = pan.condense();
                    if (condensate.length >= minimumOutputBinLength) {
                        writer.writeRecord(condensate);
                    }
                    pan.startOver(record);
                }
            }
            final IntervalsListRecord condensate = pan.condense();
            if (condensate.length >= minimumOutputBinLength) {
                writer.writeRecord(condensate);
            }
            writer.writeHeaderIfApplies();
            return "SUCCESS";
        } catch (final IOException ex) {
            throw new GATKException(ex.getMessage(), ex);
        }
    }


    private class IntervalsListReader extends TableReader<IntervalsListRecord> {
        private IntervalsListReader(final String input) throws IOException {
            super(input, new InputStreamReader(BucketUtils.openFile(input)));
        }


        private List<String> headerLines;

        @Override
        public boolean isColumnNameLine(final String[] line) {
            return super.isColumnNameLine(line)
                    || (line.length >= 3
                        && (!NumberUtils.isNumber(line[1]) || !NumberUtils.isNumber(line[2])));
        }

        @Override
        public boolean isCommentLine(final String[] line) {
            return line.length > 0 &&
                    (line[0].startsWith("#") || line[0].startsWith("@"));
        }

        @Override
        public void processCommentLine(final String line, final long number) {
            if (headerLines == null) {
                headerLines = new ArrayList<>();
            }
            headerLines.add(line);
        }

        public List<String> headerLines() {
            if (headerLines == null) {
                return Collections.emptyList();
            } else {
                return Collections.unmodifiableList(headerLines);
            }
        }

        public boolean columnNameLineIsPresent() {
            return this.columnNamesArePresent;
        }

        @Override
        public TableColumnCollection processColumns(final TableColumnCollection inputColumns) {
            if (inputColumns.columnCount() < 3) {
                throw formatException("missing compulsory coordinates columns");
            }
            return inputColumns;
        }

        @Override
        protected TableColumnCollection defaultColumns(final String[] values) {
            if (values.length < 3) {
                throw formatException("no enough elements");
            }
            final String[] columnNames = new String[values.length];
            final IntervalsListColumns[] standardColumns = IntervalsListColumns.values();
            int i;
            for (i = 0; i < columnNames.length && i < standardColumns.length; i++) {
                columnNames[i] = standardColumns[i].name();
            }
            for (int j = i; j < columnNames.length; j++) {
                columnNames[i] = String.format("_unknown_%4d", j - i);
            }
            return new TableColumnCollection(columnNames);
        }

        public IntervalsListRecord createRecord(final DataLine dataLine) {
            return new IntervalsListRecord(dataLine.get(0), dataLine.getInt(1), dataLine.getInt(2), dataLine.toArray(3));
        }
    }

    private static class IntervalsListWriter extends TableWriter<IntervalsListRecord> {

        private List<String> headerLines = new ArrayList<>();
        private boolean mustOutputColumnNameLine = false;

        private IntervalsListWriter(final String output, final TableColumnCollection columns) throws IOException {
            super(composeWriter(output), columns);
        }

        public void addHeaderLines(final Collection<String> lines) {
            headerLines.addAll(lines);
        }

        public void mustOutputColumnNameLine(final boolean value) {
            mustOutputColumnNameLine = value;
        }

        private static OutputStreamWriter composeWriter(final String output) {
            if (IOUtil.hasBlockCompressedExtension(output)) {
                return new OutputStreamWriter(new BlockCompressedOutputStream(BucketUtils.createFile(output), (File) null));
            } else {
                return new OutputStreamWriter(BucketUtils.createFile(output));
            }
        }


        @Override
        public void writeHeader() throws IOException {
            for (final String headerLine : headerLines) {
                writeLineBypassingCSVWriter(headerLine);
            }
            if (this.mustOutputColumnNameLine) {
                writeLineBypassingCSVWriter(String.join(TableUtils.COLUMN_SEPARATOR_STRING, columns.names()));
            }
        }

        public void composeLine(final IntervalsListRecord record, final DataLine line) {
            line.append(record.contig)
                    .append(record.start)
                    .append(record.end)
                    .append(record.others);
        }
    }

    private static class IntervalsListRecord {
        public final String contig;
        public final int start;
        public final int end;
        public String[] others;
        public final int length;
        private IntervalsListRecord(final String contig, final int start, final int end, String[] others) {
            this.contig = contig;
            this.start = start;
            this.end = end;
            this.others = others;
            this.length = end - start + 1;
        }

    }

    private static class CondensationPan {

        private static final String STRAND_COLUMN_NAME = "strand";
        private final TableColumnCollection columns;
        private final int condensationFactor;
        private final List<IntervalsListRecord> records;

        private CondensationPan(final int condensationFactor, final TableColumnCollection columns,final IntervalsListRecord record) {
            this.condensationFactor = condensationFactor;
            records = new ArrayList<>(Math.min(condensationFactor, 1000));
            this.columns = columns;
            records.add(record);
        }

        public static CondensationPan start(final int condensationFactor, final TableColumnCollection columns, final IntervalsListRecord record) {
            ParamUtils.isPositive(condensationFactor, "the condensation factor must be greater than 0");
            Utils.nonNull(record);
            return new CondensationPan(condensationFactor, columns, record);
        }

        public void startOver(final IntervalsListRecord record) {
            records.clear();
            records.add(record);
        }

        public boolean add(final IntervalsListRecord record) {
            if (canAdd(record)) {
                records.add(record);
                return true;
            } else {
                return false;
            }
        }

        public boolean canAdd(final IntervalsListRecord record) {
            if (records.size() >= condensationFactor) {
                return false;
            } else if (!records.get(0).contig.equals(record.contig)) {
                return false;
            } else if (record.start < records.get(records.size() - 1).start) {
                throw new UserException.BadInput("input intervals out of order: " + record.contig + ":" + record.start);
            } else {
                return record.start - records.get(records.size() - 1).end <= 1;
            }
        }

        public IntervalsListRecord condense() {
            final IntervalsListRecord first = records.get(0);
            final String contig = first.contig;
            final int start = first.start;
            final String[] summary = new String[first.others.length];

            // we skip the first 3 columns that must be contig start end.
            for (int c = 3; c < columns.columnCount(); c++) {
                final String columnName = columns.nameAt(c);
                if (columnName.equals(IntervalsListColumns.STRAND.name())) {
                    summary[c - 3] = summarizeStrand(c - 3);
                } else if (columnName.equals(IntervalsListColumns.TARGET.name())) {
                    summary[c - 3] = summarizeTarget(c - 3);
                } else {
                    summary[c - 3] = summarizeUnknownColumn(c - 3);
                }
            }
            return new IntervalsListRecord(contig, start, records.stream().mapToInt(r -> r.end).max().getAsInt() ,summary);

        }

        /**
         * If all records have the same strand we output that, otherwise we output
         * ".". If any strand is
         * @param index on the other column where the strand are stored.
         * @return never {@code null}.
         */
        private String summarizeStrand(final int index) {
            final String candidate = records.get(0).others[index];
            if (!candidate.equals(".")) {
                for (int i = 1; i < records.size(); i++) {
                    if (!records.get(i).others[index].equals(candidate)) {
                        return "."; // no strand in BED format.
                    }
                }
            }
            return candidate;
        }

        private String summarizeTarget(final int index) {
            String first = null;
            String last = null;
            boolean lastIsFirst = true;
            for (final IntervalsListRecord record : records) {
                final String name = record.others[index];
                if (name != null && !name.equals(".") && !name.trim().isEmpty()) {
                    if (first == null) {
                        first = last = name;
                        lastIsFirst = true;
                    } else {
                        last = name;
                        lastIsFirst = false;
                    }
                }
            }
            if (first == null) {
                return ".";
            } else if (lastIsFirst) {
                return first;
            } else {
                return first.trim() + "_to_"  + last.trim();
            }
        }

        private String summarizeUnknownColumn(final int index) {
            boolean seemsNumeric = true;
            boolean seemsInteger = true;
            boolean sameValue = true;
            String first = null;
            String last = null;
            double sum = 0;
                for (final IntervalsListRecord record : records) {
                    final String value = record.others[index];
                    if (value != null && !value.trim().isEmpty() && !value.equals(".")) {
                        if (seemsNumeric) {
                            if (NumberUtils.isNumber(value)) {
                                sum += Double.parseDouble(value);
                            } else {
                                seemsNumeric = false;
                            }
                        }
                        if (first == null) {
                            first = last = value;
                            sameValue = true;
                        } else {
                            last = value;
                            sameValue &= first.equals(last);
                        }
                    }
                }
            if (first == null) {
                return ".";
            } else if (seemsNumeric) {
                return Math.round(sum) == sum ? String.valueOf((long) sum) : String.valueOf(sum);
            } else if (sameValue) {
                return first;
            } else {
                return first + " ... " + last;
            }
        }
    }
}
