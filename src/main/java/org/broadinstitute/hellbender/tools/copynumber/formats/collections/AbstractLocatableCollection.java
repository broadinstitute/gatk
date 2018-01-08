package org.broadinstitute.hellbender.tools.copynumber.formats.collections;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.OverlapDetector;
import org.broadinstitute.hellbender.tools.copynumber.arguments.CopyNumberArgumentValidationUtils;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.LocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.Metadata;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;

import java.io.File;
import java.util.Comparator;
import java.util.List;
import java.util.function.BiConsumer;
import java.util.function.Function;
import java.util.stream.Collectors;

/**
 * Represents a sequence dictionary, an immutable, coordinate-sorted (with no overlaps allowed) collection of records
 * that extend {@link Locatable} (although contigs are assumed to be non-null when writing to file),
 * a set of mandatory column headers given by a {@link TableColumnCollection},
 * and lambdas for reading and writing records. Records are sorted using the sequence dictionary.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public abstract class AbstractLocatableCollection<METADATA extends LocatableMetadata, RECORD extends Locatable> extends AbstractRecordCollection<METADATA, RECORD> {
    /**
     * @param metadata records are sorted using the contained {@link SAMSequenceDictionary}
     */
    AbstractLocatableCollection(final METADATA metadata,
                                final List<RECORD> records,
                                final TableColumnCollection mandatoryColumns,
                                final Function<DataLine, RECORD> recordFromDataLineDecoder,
                                final BiConsumer<RECORD, DataLine> recordToDataLineEncoder) {
        super(metadata, sortRecords(records, metadata.getSequenceDictionary()), mandatoryColumns, recordFromDataLineDecoder, recordToDataLineEncoder);
        CopyNumberArgumentValidationUtils.validateIntervals(getRecords(), metadata.getSequenceDictionary());
    }

    /**
     * @throws IllegalArgumentException if records are not sorted according to the {@link SAMSequenceDictionary} contained in the input file
     */
    AbstractLocatableCollection(final File inputFile,
                                final TableColumnCollection mandatoryColumns,
                                final Function<DataLine, RECORD> recordFromDataLineDecoder,
                                final BiConsumer<RECORD, DataLine> recordToDataLineEncoder) {
        super(inputFile, mandatoryColumns, recordFromDataLineDecoder, recordToDataLineEncoder);
        CopyNumberArgumentValidationUtils.validateIntervals(getRecords(), getMetadata().getSequenceDictionary());
    }

    private static <T extends Locatable> List<T> sortRecords(final List<T> records,
                                                             final SAMSequenceDictionary sequenceDictionary) {
        Utils.nonNull(records);
        Utils.nonNull(sequenceDictionary);
        return records.stream()
                .sorted(IntervalUtils.getDictionaryOrderComparator(sequenceDictionary))
                .collect(Collectors.toList());
    }

    /**
     * @return  a new modifiable list of {@link SimpleInterval}s corresponding to the {@link Locatable}s
     *          for each record contained in the collection
     */
    public List<SimpleInterval> getIntervals() {
        return getRecords().stream()
                .map(r -> new SimpleInterval(r.getContig(), r.getStart(), r.getEnd()))
                .collect(Collectors.toList());
    }

    public OverlapDetector<RECORD> getOverlapDetector() {
        return OverlapDetector.create(getRecords());
    }

    public Comparator<Locatable> getComparator() {
        return IntervalUtils.getDictionaryOrderComparator(getMetadata().getSequenceDictionary());
    }

    @Override
    protected Metadata.Type getMetadataType() {
        return Metadata.Type.LOCATABLE;
    }
}