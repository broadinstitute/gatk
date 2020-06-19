package org.broadinstitute.hellbender.tools.variantdb;

import com.google.cloud.bigquery.FieldValue;
import com.google.cloud.bigquery.FieldValueList;
import com.google.cloud.bigquery.TableResult;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.utils.bigquery.BigQueryUtils;
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
        String fqSampleTableName = sampleTableRef.getFQTableName();
        return new HashSet<String>(getSampleIdMap(fqSampleTableName, printDebugInformation).values());
    }

    public static Map<Integer, String> getSampleIdMap(String fqSampleTableName, boolean printDebugInformation) {
    
        Map<Integer, String> results = new HashMap<>();

        // Get the query string:
        final String sampleListQueryString = 
            "SELECT " + SchemaUtils.SAMPLE_ID_FIELD_NAME + ", " + SchemaUtils.SAMPLE_NAME_FIELD_NAME + 
            " FROM `" + fqSampleTableName + "`";
        

        // Execute the query:
        final TableResult result = BigQueryUtils.executeQuery(sampleListQueryString);

        // Show our pretty results:
        if (printDebugInformation) {
            logger.info("Sample names returned:");
            final String prettyQueryResults = BigQueryUtils.getResultDataPrettyString(result);
            logger.info("\n" + prettyQueryResults);
        }

        // Add our samples to our map:
        for (final FieldValueList row : result.iterateAll()) {
            results.put((int) row.get(0).getLongValue(), row.get(1).getStringValue());
        }

        return results;
    }

    private static String getOptionalString(FieldValue v) {
        return (v == null || v.isNull()) ? null : v.getStringValue();
    }

    public static Map<Long, ProbeInfo> getProbeIdMap(String fqProbeTableName, boolean printDebugInformation) {
    
        Map<Long, ProbeInfo> results = new HashMap<>();

        // Get the query string:
        final String sampleListQueryString = 
            "SELECT probeId, Name, Chr, Position, Ref, AlleleA, AlleleB" + 
            " FROM `" + fqProbeTableName + "`";
        

        // Execute the query:
        final TableResult result = BigQueryUtils.executeQuery(sampleListQueryString);

        System.out.println("Beginning probe retrieval...");
        for (final FieldValueList row : result.iterateAll()) {
            ProbeInfo p = new ProbeInfo(row.get(0).getLongValue(),
                                        getOptionalString(row.get(1)), // name
                                        row.get(2).getStringValue(), // contig
                                        row.get(3).getLongValue(),   // position
                                        getOptionalString(row.get(4)), // ref
                                        getOptionalString(row.get(5)), // alleleA
                                        getOptionalString(row.get(6)));// alleleB

            results.put(p.probeId, p);
            
        }
        System.out.println("Done probe retrieval...");

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

    final static Comparator<GenericRecord> COMPRESSED_PROBE_ID_COMPARATOR = new Comparator<GenericRecord>() {
        @Override
        public int compare( GenericRecord o1, GenericRecord o2 ) {
            final long firstProbeId = RawArrayData.decode((Long) o1.get(SchemaUtils.RAW_ARRAY_DATA_FIELD_NAME)).probeId;
            final long secondProbeId = RawArrayData.decode((Long) o2.get(SchemaUtils.RAW_ARRAY_DATA_FIELD_NAME)).probeId;

            return Long.compare(firstProbeId, secondProbeId);
        }
    };

    final static Comparator<GenericRecord> UNCOMPRESSED_PROBE_ID_COMPARATOR = new Comparator<GenericRecord>() {
        @Override
        public int compare( GenericRecord o1, GenericRecord o2 ) {
            final long firstProbeId = (Long) o1.get("probe_id");
            final long secondProbeId = (Long) o2.get("probe_id");
            return Long.compare(firstProbeId, secondProbeId);
        }
    };


}
