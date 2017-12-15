package org.broadinstitute.hellbender.tools.copynumber.formats.collections;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.Metadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SampleLocatableMetadata;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;

import java.io.File;
import java.util.List;
import java.util.function.BiConsumer;
import java.util.function.Function;

/**
 * Represents a sample name, a sequence dictionary,
 * an immutable, coordinate-sorted (with no overlaps allowed) collection of records
 * that extend {@link Locatable} (although contigs are assumed to be non-null when writing to file),
 * a set of mandatory column headers given by a {@link TableColumnCollection},
 * and lambdas for reading and writing records.  Records are sorted using the sequence dictionary. See the unit test
 * for a simple example of a subclass.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public abstract class AbstractSampleLocatableCollection<RECORD extends Locatable> extends AbstractLocatableCollection<SampleLocatableMetadata, RECORD> {
    /**
     * @param metadata records are sorted using the contained {@link SAMSequenceDictionary}
     */
    AbstractSampleLocatableCollection(final SampleLocatableMetadata metadata,
                                      final List<RECORD> records,
                                      final TableColumnCollection mandatoryColumns,
                                      final Function<DataLine, RECORD> recordFromDataLineDecoder,
                                      final BiConsumer<RECORD, DataLine> recordToDataLineEncoder) {
        super(metadata, records, mandatoryColumns, recordFromDataLineDecoder, recordToDataLineEncoder);
    }

    /**
     * @throws IllegalArgumentException if records are not sorted according to the {@link SAMSequenceDictionary} contained in the input file
     */
    AbstractSampleLocatableCollection(final File inputFile,
                                      final TableColumnCollection mandatoryColumns,
                                      final Function<DataLine, RECORD> recordFromDataLineDecoder,
                                      final BiConsumer<RECORD, DataLine> recordToDataLineEncoder) {
        super(inputFile, mandatoryColumns, recordFromDataLineDecoder, recordToDataLineEncoder);
    }

    @Override
    protected Metadata.Type getMetadataType() {
        return Metadata.Type.SAMPLE_LOCATABLE;
    }
}