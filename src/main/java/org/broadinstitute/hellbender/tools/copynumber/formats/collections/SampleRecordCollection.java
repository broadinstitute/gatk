package org.broadinstitute.hellbender.tools.copynumber.formats.collections;

import com.google.common.collect.ImmutableList;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SampleMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SampleNameUtils;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SimpleSampleMetadata;
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
 * Represents an immutable collection of records associated with a sample, a set of
 * mandatory column headers given by a {@link TableColumnCollection}, and lambdas for
 * reading and writing records.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public abstract class SampleRecordCollection<RECORD> implements SampleMetadata {
    private final SampleMetadata sampleMetadata;
    private final ImmutableList<RECORD> records;
    private final TableColumnCollection mandatoryColumns;
    private final Function<DataLine, RECORD> recordFromDataLineDecoder;
    private final BiConsumer<RECORD, DataLine> recordToDataLineEncoder;

    /**
     * Constructor given the sample metadata, the list of records, the mandatory column headers,
     * and the lambdas for reading and writing records.
     *
     * @param sampleMetadata                metadata associated with the sample
     * @param records                       list of records associated with the sample; may be empty
     * @param mandatoryColumns              mandatory columns required to construct collection from a TSV file; cannot be empty
     * @param recordFromDataLineDecoder     lambda for decoding a record from a {@link DataLine} when reading from a TSV file
     * @param recordToDataLineEncoder       lambda for encoding a record to a {@link DataLine} when writing to a TSV file
     */
    SampleRecordCollection(final SampleMetadata sampleMetadata,
                           final List<RECORD> records,
                           final TableColumnCollection mandatoryColumns,
                           final Function<DataLine, RECORD> recordFromDataLineDecoder,
                           final BiConsumer<RECORD, DataLine> recordToDataLineEncoder) {
        this.sampleMetadata = Utils.nonNull(sampleMetadata);
        this.records = ImmutableList.copyOf(Utils.nonNull(records));
        this.mandatoryColumns = Utils.nonNull(mandatoryColumns);
        this.recordFromDataLineDecoder = Utils.nonNull(recordFromDataLineDecoder);
        this.recordToDataLineEncoder = Utils.nonNull(recordToDataLineEncoder);
        Utils.nonEmpty(mandatoryColumns.names());
    }

    /**
     * Constructor given an input file, the mandatory column headers, and the lambdas for reading and writing records.
     * The sample metadata is read from the header of the file and the list of records is read using
     * the column headers and the appropriate lambda.
     *
     * @param inputFile                     TSV file containing sample metadata and records; the latter may be empty
     * @param mandatoryColumns              mandatory columns required to construct collection from a TSV file; cannot be empty
     * @param recordFromDataLineDecoder     lambda for decoding a record from a {@link DataLine} when reading from a TSV file
     * @param recordToDataLineEncoder       lambda for encoding a record to a {@link DataLine} when writing to a TSV file
     */
    SampleRecordCollection(final File inputFile,
                           final TableColumnCollection mandatoryColumns,
                           final Function<DataLine, RECORD> recordFromDataLineDecoder,
                           final BiConsumer<RECORD, DataLine> recordToDataLineEncoder) {
        IOUtils.canReadFile(inputFile);
        this.mandatoryColumns = Utils.nonNull(mandatoryColumns);
        this.recordFromDataLineDecoder = Utils.nonNull(recordFromDataLineDecoder);
        this.recordToDataLineEncoder = Utils.nonNull(recordToDataLineEncoder);
        Utils.nonEmpty(mandatoryColumns.names());

        try (final SampleRecordCollectionReader reader = new SampleRecordCollectionReader(inputFile)) {
            TableUtils.checkMandatoryColumns(reader.columns(), mandatoryColumns, UserException.BadInput::new);
            sampleMetadata = reader.readSampleMetadata();
            records = reader.stream().collect(Collectors.collectingAndThen(Collectors.toList(), ImmutableList::copyOf));
        } catch (final IOException | UncheckedIOException e) {
            throw new UserException.CouldNotReadInputFile(inputFile, e);
        }
    }

    public final SampleMetadata getSampleMetadata() {
        return sampleMetadata;
    }

    @Override
    public final String getSampleName() {
        return sampleMetadata.getSampleName();
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
     * Writes the sample metadata and the records to file.
     */
    public void write(final File outputFile) {
        try (final SampleRecordCollectionWriter writer = new SampleRecordCollectionWriter(outputFile, sampleMetadata)) {
            writer.writeSampleMetadata();
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

        final SampleRecordCollection<?> that = (SampleRecordCollection<?>) o;
        return sampleMetadata.equals(that.sampleMetadata) && records.equals(that.records);
    }

    @Override
    public int hashCode() {
        int result = sampleMetadata.hashCode();
        result = 31 * result + records.hashCode();
        return result;
    }

    private final class SampleRecordCollectionReader extends TableReader<RECORD> {
        private final File file;

        private SampleRecordCollectionReader(final File file) throws IOException {
            super(file);
            this.file = file;
        }

        private SampleMetadata readSampleMetadata() {
            final String sampleName = SampleNameUtils.readSampleName(file);
            return new SimpleSampleMetadata(sampleName);
        }

        @Override
        protected RECORD createRecord(final DataLine dataLine) {
            Utils.nonNull(dataLine);
            return recordFromDataLineDecoder.apply(dataLine);
        }
    }

    private final class SampleRecordCollectionWriter extends TableWriter<RECORD> {
        private final SampleMetadata sampleMetadata;

        SampleRecordCollectionWriter(final File file,
                                     final SampleMetadata sampleMetadata) throws IOException {
            super(file, mandatoryColumns);
            this.sampleMetadata = Utils.nonNull(sampleMetadata);
        }

        void writeSampleMetadata() {
            writeSampleName();
        }

        private void writeSampleName() {
            try {
                writeComment(SampleNameUtils.SAMPLE_NAME_COMMENT_TAG + sampleMetadata.getSampleName());
            } catch (final IOException e) {
                throw new UserException("Could not write sample name.");
            }
        }

        @Override
        protected void composeLine(final RECORD record, final DataLine dataLine) {
            Utils.nonNull(record);
            Utils.nonNull(dataLine);
            recordToDataLineEncoder.accept(record, dataLine);
        }
    }
}