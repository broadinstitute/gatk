package org.broadinstitute.hellbender.tools.variantdb.arrays;

import com.google.cloud.bigquery.FieldValue;
import com.google.cloud.bigquery.FieldValueList;
import com.google.cloud.bigquery.TableResult;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.tools.variantdb.SchemaUtils;
import org.broadinstitute.hellbender.tools.variantdb.arrays.tables.ProbeInfo;
import org.broadinstitute.hellbender.tools.variantdb.arrays.tables.SampleList;
import org.broadinstitute.hellbender.utils.bigquery.TableReference;
import org.broadinstitute.hellbender.utils.localsort.AvroSortingCollectionCodec;
import org.broadinstitute.hellbender.utils.localsort.SortingCollection;
import org.apache.avro.generic.GenericRecord;

import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

public class ExtractCohortBQ {
    private static final Logger logger = LogManager.getLogger(ExtractCohortBQ.class);

    public static Set<String> populateSampleNames(TableReference sampleTableRef, boolean printDebugInformation) {
        return new HashSet<>(SampleList.getSampleIdMap(sampleTableRef, printDebugInformation).values());
    }

    public static Map<String, ProbeInfo> getProbeNameMap(String fqProbeTableName, boolean printDebugInformation) {
        Map<String, ProbeInfo> results = new HashMap<>();
        for (final ProbeInfo pi : ProbeInfo.getProbeIdMapWithStorageAPI(fqProbeTableName, printDebugInformation).values()) {
            results.put(pi.name, pi);
        }
        return results;
    }

    public static SortingCollection<GenericRecord> getAvroSortingCollection(org.apache.avro.Schema schema, int localSortMaxRecordsInRam) {
        final SortingCollection.Codec<GenericRecord> sortingCollectionCodec = new AvroSortingCollectionCodec(schema);
        final Comparator<GenericRecord> sortingCollectionComparator = new Comparator<GenericRecord>() {
            @Override
            public int compare( GenericRecord o1, GenericRecord o2 ) {
                final long firstPosition = Long.parseLong(o1.get(SchemaUtils.LOCATION_FIELD_NAME).toString());
                final long secondPosition = Long.parseLong(o2.get(SchemaUtils.LOCATION_FIELD_NAME).toString());

                return Long.compare(firstPosition, secondPosition);
            }
        };
        return SortingCollection.newInstance(GenericRecord.class, sortingCollectionCodec, sortingCollectionComparator, localSortMaxRecordsInRam);
    }

    public static SortingCollection<GenericRecord> getAvroProbeIdSortingCollection(org.apache.avro.Schema schema, int localSortMaxRecordsInRam, Comparator<GenericRecord> comparator) {
        final SortingCollection.Codec<GenericRecord> sortingCollectionCodec = new AvroSortingCollectionCodec(schema);
        return SortingCollection.newInstance(GenericRecord.class, sortingCollectionCodec, comparator, localSortMaxRecordsInRam);
    }

    public final static Comparator<GenericRecord> COMPRESSED_PROBE_ID_COMPARATOR = new Comparator<GenericRecord>() {
        @Override
        public int compare( GenericRecord o1, GenericRecord o2 ) {
            final long firstProbeId = new BasicArrayData((Long) o1.get(SchemaUtils.BASIC_ARRAY_DATA_FIELD_NAME)).probeId;
            final long secondProbeId = new BasicArrayData((Long) o2.get(SchemaUtils.BASIC_ARRAY_DATA_FIELD_NAME)).probeId;

            return Long.compare(firstProbeId, secondProbeId);
        }
    };

    public final static Comparator<GenericRecord> UNCOMPRESSED_PROBE_ID_COMPARATOR = new Comparator<GenericRecord>() {
        @Override
        public int compare( GenericRecord o1, GenericRecord o2 ) {
            final long firstProbeId = (Long) o1.get("probe_id");
            final long secondProbeId = (Long) o2.get("probe_id");
            return Long.compare(firstProbeId, secondProbeId);
        }
    };


}
