package org.broadinstitute.hellbender.tools.gvs.ingest.bq;

import org.broadinstitute.hellbender.tools.gvs.common.VrsAlleleRecord;
import org.broadinstitute.hellbender.tools.gvs.ingest.VrsAlleleWriter;

import java.io.IOException;

/**
 * BigQuery writer for the {@code vrs_allele} table.
 */
public class VrsAlleleBQWriter extends AbstractBQWriter implements VrsAlleleWriter {

    public VrsAlleleBQWriter(String projectId, String datasetName, String tableName) throws IOException {
        super(projectId, datasetName, tableName);
    }

    @Override
    public void write(VrsAlleleRecord record) throws IOException {
        write(toJson(record));
    }
}
