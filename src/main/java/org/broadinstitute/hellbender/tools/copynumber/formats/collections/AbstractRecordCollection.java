package org.broadinstitute.hellbender.tools.copynumber.formats.collections;

import com.google.common.collect.ImmutableList;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMTextHeaderCodec;
import htsjdk.samtools.util.BufferedLineReader;
import htsjdk.samtools.util.LineReader;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.formats.CopyNumberFormatsUtils;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.Metadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.MetadataUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.tsv.*;

import java.io.*;
import java.util.List;
import java.util.function.BiConsumer;
import java.util.function.Function;
import java.util.stream.Collectors;

/**
 * Represents {@link METADATA} (which can be represented as a {@link SAMFileHeader}),
 * an immutable collection of records,
 * a set of mandatory column headers given by a {@link TableColumnCollection},
 * and lambdas for reading and writing records.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public abstract class AbstractRecordCollection<METADATA extends Metadata, RECORD> {
    private final METADATA metadata;
    private final ImmutableList<RECORD> records;
    private final TableColumnCollection mandatoryColumns;
    private final Function<DataLine, RECORD> recordFromDataLineDecoder;
    private final BiConsumer<RECORD, DataLine> recordToDataLineEncoder;

    /**
     * Constructor given the {@link METADATA}, the list of records, the mandatory column headers,
     * and the lambdas for reading and writing records.
     *
     * @param metadata                      {@link METADATA} (which can be represented as a {@link SAMFileHeader}
     * @param records                       list of records; may be empty
     * @param mandatoryColumns              mandatory columns required to construct collection from a TSV file; cannot be empty
     * @param recordFromDataLineDecoder     lambda for decoding a record from a {@link DataLine} when reading from a TSV file
     * @param recordToDataLineEncoder       lambda for encoding a record to a {@link DataLine} when writing to a TSV file
     */
    AbstractRecordCollection(final METADATA metadata,
                             final List<RECORD> records,
                             final TableColumnCollection mandatoryColumns,
                             final Function<DataLine, RECORD> recordFromDataLineDecoder,
                             final BiConsumer<RECORD, DataLine> recordToDataLineEncoder) {
        this.metadata = Utils.nonNull(metadata);
        this.records = ImmutableList.copyOf(Utils.nonNull(records));
        this.mandatoryColumns = Utils.nonNull(mandatoryColumns);
        this.recordFromDataLineDecoder = Utils.nonNull(recordFromDataLineDecoder);
        this.recordToDataLineEncoder = Utils.nonNull(recordToDataLineEncoder);
        Utils.nonEmpty(mandatoryColumns.names());
    }

    /**
     * Constructor given an input file, the mandatory column headers, and the lambdas for reading and writing records.
     * The list of records is read using the column headers and the appropriate lambda.
     *
     * @param inputFile                     TSV file; must contain a {@link SAMFileHeader} and mandatory column headers, but can contain no records
     * @param mandatoryColumns              mandatory columns required to construct collection from a TSV file; cannot be empty
     * @param recordFromDataLineDecoder     lambda for decoding a record from a {@link DataLine} when reading from a TSV file
     * @param recordToDataLineEncoder       lambda for encoding a record to a {@link DataLine} when writing to a TSV file
     */
    AbstractRecordCollection(final File inputFile,
                             final TableColumnCollection mandatoryColumns,
                             final Function<DataLine, RECORD> recordFromDataLineDecoder,
                             final BiConsumer<RECORD, DataLine> recordToDataLineEncoder) {
        IOUtils.canReadFile(inputFile);
        this.mandatoryColumns = Utils.nonNull(mandatoryColumns);
        this.recordFromDataLineDecoder = Utils.nonNull(recordFromDataLineDecoder);
        this.recordToDataLineEncoder = Utils.nonNull(recordToDataLineEncoder);
        Utils.nonEmpty(mandatoryColumns.names());

        try (final RecordCollectionReader reader = new RecordCollectionReader(inputFile)) {
            metadata = MetadataUtils.fromHeader(reader.getHeader(), getMetadataType());
            TableUtils.checkMandatoryColumns(reader.columns(), mandatoryColumns, UserException.BadInput::new);
            records = reader.stream().collect(Collectors.collectingAndThen(Collectors.toList(), ImmutableList::copyOf));
        } catch (final IOException | UncheckedIOException e) {
            throw new UserException.CouldNotReadInputFile(inputFile, e);
        }
    }

    public final int size() {
        return records.size();
    }

    /**
     * Subclasses should add an enum to {@link Metadata.Type}, a corresponding switch case statement
     * to {@link MetadataUtils#fromHeader(SAMFileHeader, Metadata.Type)}, and implement this method accordingly.
     */
    abstract Metadata.Type getMetadataType();

    public METADATA getMetadata() {
        return metadata;
    }

    /**
     * @return  an immutable view of the records contained in the collection
     */
    public final List<RECORD> getRecords() {
        return records;
    }

    /**
     * Writes the records to file.
     */
    public void write(final File outputFile) {
        try (final FileWriter writer = new FileWriter(outputFile)) {
            writer.write(metadata.toHeader().getSAMString());
        } catch (final IOException e) {
            throw new UserException.CouldNotCreateOutputFile(outputFile, e);
        }
        try (final RecordWriter recordWriter = new RecordWriter(new FileWriter(outputFile, true))) {
            recordWriter.writeAllRecords(records);
        } catch (final IOException e) {
            throw new UserException.CouldNotCreateOutputFile(outputFile, e);
        }
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) {
            return true;
        }
        if (o == null || getClass() != o.getClass()) {
            return false;
        }

        final AbstractRecordCollection<?, ?> that = (AbstractRecordCollection<?, ?>) o;

        return metadata.equals(that.metadata) &&
                records.equals(that.records) &&
                mandatoryColumns.equals(that.mandatoryColumns) &&
                recordFromDataLineDecoder.equals(that.recordFromDataLineDecoder) &&
                recordToDataLineEncoder.equals(that.recordToDataLineEncoder);
    }

    @Override
    public int hashCode() {
        int result = metadata.hashCode();
        result = 31 * result + records.hashCode();
        result = 31 * result + mandatoryColumns.hashCode();
        result = 31 * result + recordFromDataLineDecoder.hashCode();
        result = 31 * result + recordToDataLineEncoder.hashCode();
        return result;
    }

    @Override
    public String toString() {
        return "AbstractRecordCollection{" +
                "metadata=" + metadata +
                ", records=" + records +
                '}';
    }

    static String formatDouble(final double value) {
        return CopyNumberFormatsUtils.formatDouble(value);
    }

    final class RecordCollectionReader extends TableReader<RECORD> {
        private static final String COMMENT_PREFIX = CopyNumberFormatsUtils.COMMENT_PREFIX;   //SAMTextHeaderCodec.HEADER_LINE_START; we need TableReader to treat SAM header as comment lines
        private final File file;

        RecordCollectionReader(final File file) throws IOException {
            super(file);
            this.file = file;
        }

        @Override
        protected RECORD createRecord(final DataLine dataLine) {
            Utils.nonNull(dataLine);
            return recordFromDataLineDecoder.apply(dataLine);
        }

        private SAMFileHeader getHeader() throws FileNotFoundException {
            final LineReader lineReader = new BufferedLineReader(new FileInputStream(file));
            return new SAMTextHeaderCodec().decode(lineReader, getSource());
        }

        @Override
        protected boolean isCommentLine(final String[] line) {
            return line.length > 0 && line[0].startsWith(COMMENT_PREFIX);
        }
    }

    final class RecordWriter extends TableWriter<RECORD> {
        RecordWriter(final Writer writer) throws IOException {
            super(writer, mandatoryColumns);
        }

        @Override
        protected void composeLine(final RECORD record, final DataLine dataLine) {
            Utils.nonNull(record);
            Utils.nonNull(dataLine);
            recordToDataLineEncoder.accept(record, dataLine);
        }
    }
}