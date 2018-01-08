package org.broadinstitute.hellbender.tools.copynumber.formats.collections;

import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.Metadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SampleMetadata;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;

import java.io.File;
import java.util.List;
import java.util.function.BiConsumer;
import java.util.function.Function;

/**
 * Represents a sample name,
 * an immutable collection of records,
 * a set of mandatory column headers given by a {@link TableColumnCollection},
 * and lambdas for reading and writing records.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public abstract class AbstractSampleRecordCollection<RECORD> extends AbstractRecordCollection<SampleMetadata, RECORD> {
    AbstractSampleRecordCollection(final SampleMetadata metadata,
                                   final List<RECORD> records,
                                   final TableColumnCollection mandatoryColumns,
                                   final Function<DataLine, RECORD> recordFromDataLineDecoder,
                                   final BiConsumer<RECORD, DataLine> recordToDataLineEncoder) {
        super(metadata, records, mandatoryColumns, recordFromDataLineDecoder, recordToDataLineEncoder);
    }

    AbstractSampleRecordCollection(final File inputFile,
                                   final TableColumnCollection mandatoryColumns,
                                   final Function<DataLine, RECORD> recordFromDataLineDecoder,
                                   final BiConsumer<RECORD, DataLine> recordToDataLineEncoder) {
        super(inputFile, mandatoryColumns, recordFromDataLineDecoder, recordToDataLineEncoder);
    }

    @Override
    protected Metadata.Type getMetadataType() {
        return Metadata.Type.SAMPLE;
    }
}