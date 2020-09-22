package org.broadinstitute.hellbender.tools.variantdb.arrays.tables;

import org.apache.avro.generic.GenericRecord;
import org.broadinstitute.hellbender.utils.bigquery.StorageAPIAvroReader;
import org.broadinstitute.hellbender.utils.bigquery.TableReference;

import java.util.HashMap;
import java.util.Map;

public class ProbeQcMetrics {
    public final long probeId;
    public final Double excess_het;
    public final Double call_rate;
    public final Boolean invariant;

    public ProbeQcMetrics(final long probeId, final Double excess_het, final Double call_rate, final Boolean invariant) {
        this.probeId = probeId;
        this.excess_het = excess_het;
        this.call_rate = call_rate;
        this.invariant = invariant;
    }

    public static Map<Long, ProbeQcMetrics> getProbeQcMetricsWithStorageAPI(String fqProbeTableName) {
        Map<Long, ProbeQcMetrics> results = new HashMap<>();

        TableReference tableRef = new TableReference(fqProbeTableName, ProbeQcMetricsSchema.PROBE_QC_METRIC_FIELDS);

        System.out.println("Beginning probe QC metrics retrieval...");
        long start = System.currentTimeMillis();

        try (final StorageAPIAvroReader reader = new StorageAPIAvroReader(tableRef)) {
            for ( final GenericRecord row : reader ) {                
                ProbeQcMetrics p = new ProbeQcMetrics(
                    (Long) row.get(ProbeQcMetricsSchema.PROBE_ID),
                    getOptionalDouble(row, ProbeQcMetricsSchema.EXCESS_HET),
                    getOptionalDouble(row, ProbeQcMetricsSchema.CALL_RATE),
                    getOptionalBoolean(row, ProbeQcMetricsSchema.INVARIANT)
                );
                results.put(p.probeId, p);
            }
        }

        long elapsed = System.currentTimeMillis() - start;
        System.out.println("Completed probe qc metrics retrieval... (" + elapsed + " ms)");

        return results;
    }

    private static Double getOptionalDouble(GenericRecord rec, String fieldName) {
        Object o = rec.get(fieldName);
        return ( o == null ? null : ((Double) o) );
    }

    private static Boolean getOptionalBoolean(GenericRecord rec, String fieldName) {
        Object o = rec.get(fieldName);
        return ( o == null ? null : ((Boolean) o) );
    }

    @Override
    public String toString() {
        return "ProbeQcMetric [probeId=" + probeId + ", excess_het=" + excess_het + ", call_rate="
                + call_rate + ", invariant=" + invariant + "]";
    }


    

}
