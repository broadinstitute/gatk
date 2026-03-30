package org.broadinstitute.hellbender.tools.gvs.ingest.tsv;

import org.broadinstitute.hellbender.tools.gvs.common.IngestConstants;
import org.broadinstitute.hellbender.tools.gvs.common.VrsAlleleRecord;
import org.broadinstitute.hellbender.tools.gvs.ingest.VrsAlleleWriter;
import org.broadinstitute.hellbender.utils.tsv.SimpleXSVWriter;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.ArrayList;
import java.util.List;

/**
 * TSV writer for the {@code vrs_allele} table.
 *
 * <p>Column order matches the BigQuery schema defined in {@code GvsCreateTables.wdl}:
 * <ol>
 *   <li>vrs_allele_id</li>
 *   <li>vrs_location_id</li>
 *   <li>refget_accession</li>
 *   <li>ref_genome_coord_start</li>
 *   <li>ref_genome_coord_end</li>
 *   <li>state_type</li>
 *   <li>state_sequence (nullable)</li>
 *   <li>state_length (nullable)</li>
 *   <li>state_repeat_subunit_length (nullable)</li>
 * </ol>
 */
public class VrsAlleleTsvWriter extends SimpleXSVWriter implements VrsAlleleWriter {

    private static final List<String> HEADERS = Arrays.asList(
        "vrs_allele_id",
        "vrs_location_id",
        "refget_accession",
        "ref_genome_coord_start",
        "ref_genome_coord_end",
        "state_type",
        "state_sequence",
        "state_length",
        "state_repeat_subunit_length"
    );

    public VrsAlleleTsvWriter(File file) throws IOException {
        super(file.toPath(), IngestConstants.SEPARATOR);
        this.setHeaderLine(HEADERS);
    }

    @Override
    public void write(VrsAlleleRecord record) throws IOException {
        List<String> row = new ArrayList<>();
        row.add(record.vrsAlleleId);
        row.add(record.vrsLocationId);
        row.add(record.refgetAccession);
        row.add(String.valueOf(record.start));
        row.add(String.valueOf(record.end));
        row.add(record.stateType);
        row.add(record.stateSequence  != null ? record.stateSequence               : "");
        row.add(record.stateLength    != null ? String.valueOf(record.stateLength)  : "");
        row.add(record.stateRepeatSubunitLength != null
                ? String.valueOf(record.stateRepeatSubunitLength) : "");
        this.getNewLineBuilder().setRow(row).write();
    }
}
