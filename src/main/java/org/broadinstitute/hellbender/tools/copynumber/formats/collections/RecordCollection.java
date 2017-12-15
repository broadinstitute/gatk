package org.broadinstitute.hellbender.tools.copynumber.formats.collections;

import com.google.common.collect.ImmutableList;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.tsv.*;

import java.io.File;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.util.List;
import java.util.function.BiConsumer;
import java.util.function.Function;
import java.util.stream.Collectors;

/**
 * Represents an immutable collection of records, a set of
 * mandatory column headers given by a {@link TableColumnCollection}, and lambdas for
 * reading and writing records.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public abstract class RecordCollection<RECORD> {
    private final ImmutableList<RECORD> records;
    private final TableColumnCollection mandatoryColumns;
    private final Function<DataLine, RECORD> recordFromDataLineDecoder;
    private final BiConsumer<RECORD, DataLine> recordToDataLineEncoder;

    /**
     * Constructor given the list of records, the mandatory column headers,
     * and the lambdas for reading and writing records.
     *
     * @param records                       list of records; may be empty
     * @param mandatoryColumns              mandatory columns required to construct collection from a TSV file; cannot be empty
     * @param recordFromDataLineDecoder     lambda for decoding a record from a {@link DataLine} when reading from a TSV file
     * @param recordToDataLineEncoder       lambda for encoding a record to a {@link DataLine} when writing to a TSV file
     */
    RecordCollection(final List<RECORD> records,
                     final TableColumnCollection mandatoryColumns,
                     final Function<DataLine, RECORD> recordFromDataLineDecoder,
                     final BiConsumer<RECORD, DataLine> recordToDataLineEncoder) {
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
     * @param inputFile                     TSV file containing records; may be empty
     * @param mandatoryColumns              mandatory columns required to construct collection from a TSV file; cannot be empty
     * @param recordFromDataLineDecoder     lambda for decoding a record from a {@link DataLine} when reading from a TSV file
     * @param recordToDataLineEncoder       lambda for encoding a record to a {@link DataLine} when writing to a TSV file
     */
    RecordCollection(final File inputFile,
                     final TableColumnCollection mandatoryColumns,
                     final Function<DataLine, RECORD> recordFromDataLineDecoder,
                     final BiConsumer<RECORD, DataLine> recordToDataLineEncoder) {
        IOUtils.canReadFile(inputFile);
        this.mandatoryColumns = Utils.nonNull(mandatoryColumns);
        this.recordFromDataLineDecoder = Utils.nonNull(recordFromDataLineDecoder);
        this.recordToDataLineEncoder = Utils.nonNull(recordToDataLineEncoder);
        Utils.nonEmpty(mandatoryColumns.names());

        try (final RecordCollectionReader reader = new RecordCollectionReader(inputFile)) {
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
     * @return  an immutable view of the records contained in the collection
     */
    public final List<RECORD> getRecords() {
        return records;
    }

    /**
     * Writes the records to file.
     */
    public void write(final File outputFile) {
        try (final RecordCollectionWriter writer = new RecordCollectionWriter(outputFile)) {
            writer.writeAllRecords(records);
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

        final RecordCollection<?> that = (RecordCollection<?>) o;
        return records.equals(that.records);
    }

    @Override
    public int hashCode() {
        return records.hashCode();
    }

    private final class RecordCollectionReader extends TableReader<RECORD> {
        private final File file;

        private RecordCollectionReader(final File file) throws IOException {
            super(file);
            this.file = file;
        }

        @Override
        protected RECORD createRecord(final DataLine dataLine) {
            Utils.nonNull(dataLine);
            return recordFromDataLineDecoder.apply(dataLine);
        }
    }

    private final class RecordCollectionWriter extends TableWriter<RECORD> {
        RecordCollectionWriter(final File file) throws IOException {
            super(file, mandatoryColumns);
        }

        @Override
        protected void composeLine(final RECORD record, final DataLine dataLine) {
            Utils.nonNull(record);
            Utils.nonNull(dataLine);
            recordToDataLineEncoder.accept(record, dataLine);
        }
    }
}