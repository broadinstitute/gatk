package org.broadinstitute.hellbender.tools.copynumber;

import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.IOUtil;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.SimpleCountCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.SimpleCount;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableReader;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;

import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.Writer;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

@CommandLineProgramProperties(summary = "condenses read-counts in .tsv format into a new file also in .tsv format",
        oneLineSummary = "condenses read-counts in .tsv format into a new file also in .tsv formas" ,
        programGroup = CopyNumberProgramGroup.class)
public class CondenseReadCounts extends CommandLineProgram {

    public static final String FACTOR_LONG_NAME = "factor";
    public static final String FACTOR_SHORT_NAME = "f";
    public static final String EXPECTED_BIN_LENGTH_LONG_NAME = "out-bin-length";

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
            doc = "Indicates how many input read count record must be condense into each output record",
            fullName = FACTOR_LONG_NAME,
            shortName = FACTOR_SHORT_NAME,
            minValue = 1
    )
    private int factor = 0;

    @Argument(
            doc = "expected bin-length after condensation, those output bin that do not match this length will be excluded. " +
                    "The default 0 means that such constrain is not be applied",
            fullName = EXPECTED_BIN_LENGTH_LONG_NAME,
            minValue = 0,
            optional = true
    )
    private int expectedBinLength = 0;

    @Override
    public String doWork() {

        try (final SimpleCountReader reader = new SimpleCountReader(input);
             final SimpleCountWriter writer = new SimpleCountWriter(output)) {
            String contig = null;
            int firstPosition = -1;
            int lastPosition = -1;
            int currentReadCount = 0;
            int currentRecordCount = 0;
            SimpleCount record;
            int skippedRecords = 0;
            double skippedReadCounts = 0;
            double outputCounts = 0;
            int condensedRecords = 0;
            while ((record = reader.readRecord()) != null) {
                if (record.getContig().equals(contig) && record.getStart() == lastPosition + 1) {
                    currentReadCount += record.getCount();
                    lastPosition = record.getEnd();
                    if (++currentRecordCount == factor) {
                        if (expectedBinLength > 0 && (lastPosition - firstPosition + 1) == expectedBinLength) {
                            if (condensedRecords == 0) {
                                writer.writeHeaderComments(reader.headerLines());
                            }
                            writer.writeRecord(composeRecord(contig, firstPosition, lastPosition, currentReadCount));
                            condensedRecords += currentRecordCount;
                            outputCounts += currentReadCount;
                        } else {
                            skippedReadCounts += currentReadCount;
                            skippedRecords += currentRecordCount;
                        }
                        contig = null;
                        currentRecordCount = 0;
                        currentReadCount = 0;
                        firstPosition = lastPosition = -1;
                    }
                } else {
                    skippedRecords += currentRecordCount;
                    skippedReadCounts += currentReadCount;
                    contig = record.getContig();
                    firstPosition = record.getStart();
                    lastPosition = record.getEnd();
                    currentReadCount = record.getCount();
                    currentRecordCount = 1;
                }

            }
            skippedRecords += currentRecordCount;
            skippedReadCounts += currentReadCount;
            logger.info("Number of output records: " + condensedRecords / factor);
            logger.info("Number of condensed input records: " + condensedRecords);
            logger.info("Number of skipped input records: " + skippedRecords);
            logger.info("Total count condensed: " + outputCounts);
            logger.info("Total count skipped: " + skippedReadCounts);
            logger.info(String.format("%% of lossed records and counts are: %.2g and %.2g",
                    100.00 *(skippedRecords/ (double)(skippedRecords + condensedRecords)),
                    100.00 *(skippedReadCounts/ (skippedReadCounts + outputCounts))));

            return "SUCCESS";
        } catch (final IOException ex) {
            throw new GATKException(ex.getMessage(), ex);
        }
    }

    private SimpleCount composeRecord(final String contig, final int start, final int end, final int count) {
        return new SimpleCount(new SimpleInterval(contig, start, end), count);
    }

    private class SimpleCountReader extends TableReader<SimpleCount> {
        private SimpleCountReader(final String input) throws IOException {
            super(input, new InputStreamReader(BucketUtils.openFile(input)));
        }

        private List<String> headerLines;
        private boolean foundHeader = false;

        @Override
        public boolean isCommentLine(String[] line) {
            return line.length > 0 &&
                    (line[0].startsWith("#") || line[0].startsWith("@"));
        }

        @Override
        public void processCommentLine(final String line, final long number) {
            if (!foundHeader) {
                if (headerLines == null) {
                    headerLines = new ArrayList<>();
                }
                headerLines.add(line);
            }
        }

        public List<String> headerLines() {
            return headerLines == null ? Collections.emptyList() : Collections.unmodifiableList(headerLines);
        }

        @Override
        public void processColumns(final TableColumnCollection inputColumns) {
            if (!inputColumns.containsAll(SimpleCountCollection.SimpleCountTableColumn.COLUMNS.names())) {
                throw formatException("missing columns");
            }
            foundHeader = true;
        }

        public SimpleCount createRecord(final DataLine dataLine) {
            return SimpleCountCollection.SIMPLE_COUNT_RECORD_FROM_DATA_LINE_DECODER.apply(dataLine);
        }
    }

    private static class SimpleCountWriter extends TableWriter<SimpleCount> {

        private SimpleCountWriter(final String output) throws IOException {
            super(composeWriter(output), SimpleCountCollection.SimpleCountTableColumn.COLUMNS);
        }

        private static OutputStreamWriter composeWriter(final String output) {
            if (IOUtil.hasBlockCompressedExtension(output)) {
                return new OutputStreamWriter(new BlockCompressedOutputStream(BucketUtils.createFile(output), (File) null));
            } else {
                return new OutputStreamWriter(BucketUtils.createFile(output));
            }
        }

        public void writeHeaderComments(final List<String> comments) throws IOException {
            for (final String line : comments) {
                writeLineBypassingCSVWriter(line);
            }
        }

        public void composeLine(final SimpleCount record, final DataLine line) {
            SimpleCountCollection.SIMPLE_COUNT_RECORD_TO_DATA_LINE_ENCODER.accept(record, line);
        }
    }
}
