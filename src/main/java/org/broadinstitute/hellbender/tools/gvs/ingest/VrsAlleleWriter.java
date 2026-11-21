package org.broadinstitute.hellbender.tools.gvs.ingest;

import org.broadinstitute.hellbender.tools.gvs.common.VrsAlleleRecord;
import org.json.JSONObject;

import java.io.Closeable;
import java.io.IOException;

/**
 * Writer interface for the {@code vrs_allele_scratch} table.
 *
 * <p>One row is written per normalized alternate allele per sample, carrying the sample id and the
 * full VRS record (allele ID, location ID, coordinates, refget accession, and state fields).
 * The same {@code vrs_allele_id} recurs across samples; {@code ProcessVrsAlleles} merges the distinct
 * rows into the deduplicated {@code vrs_allele} table.
 */
public interface VrsAlleleWriter extends Closeable {

    void write(Long sampleId, VrsAlleleRecord record) throws IOException;

    default void commitData() {
        // no-op by default
    }

    /**
     * Serialize a {@link VrsAlleleRecord} to a {@link JSONObject} matching the
     * {@code vrs_allele_scratch} table schema.
     */
    default JSONObject toJson(Long sampleId, VrsAlleleRecord record) {
        JSONObject obj = new JSONObject();
        obj.put("sample_id",                  sampleId);
        obj.put("vrs_allele_id",              record.vrsAlleleId);
        obj.put("vrs_location_id",            record.vrsLocationId);
        obj.put("refget_accession",           record.refgetAccession);
        obj.put("ref_genome_coord_start",     record.start);
        obj.put("ref_genome_coord_end",       record.end);
        obj.put("state_type",                 record.stateType);
        // nullable fields — only include when non-null
        if (record.stateSequence != null) {
            obj.put("state_sequence", record.stateSequence);
        }
        if (record.stateLength != null) {
            obj.put("state_length", record.stateLength);
        }
        if (record.stateRepeatSubunitLength != null) {
            obj.put("state_repeat_subunit_length", record.stateRepeatSubunitLength);
        }
        return obj;
    }
}
