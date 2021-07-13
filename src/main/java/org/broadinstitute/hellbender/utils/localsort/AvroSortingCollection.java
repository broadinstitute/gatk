package org.broadinstitute.hellbender.utils.localsort;

import org.apache.avro.generic.GenericRecord;

import java.util.Comparator;

public class AvroSortingCollection {

    public static SortingCollection<GenericRecord> getAvroSortingCollection(org.apache.avro.Schema schema, int localSortMaxRecordsInRam, Comparator<GenericRecord> comparator) {
        final SortingCollection.Codec<GenericRecord> sortingCollectionCodec = new AvroSortingCollectionCodec(schema);
        return SortingCollection.newInstance(GenericRecord.class, sortingCollectionCodec, comparator, localSortMaxRecordsInRam);
    }

}
