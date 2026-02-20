package org.broadinstitute.hellbender.tools.gvs.ingest.bq;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.tools.gvs.ingest.VetWriter;
import org.json.JSONObject;

import java.io.IOException;

public class VetBQWriter extends AbstractBQWriter implements VetWriter {

    private final boolean skipLoadingVqsrFields;
    private final boolean forceLoadingFromNonAlleleSpecific;

    public VetBQWriter(String projectId, String datasetName, String tableName, boolean skipLoadingVqsrFields, boolean forceLoadingFromNonAlleleSpecific) throws IOException {
        super(projectId, datasetName, tableName);
        this.skipLoadingVqsrFields = skipLoadingVqsrFields;
        this.forceLoadingFromNonAlleleSpecific = forceLoadingFromNonAlleleSpecific;
    }

    @Override
    public void write(long location, VariantContext variant, long sampleId) throws IOException {
        JSONObject jsonObject = createJson(location, variant, sampleId, skipLoadingVqsrFields, forceLoadingFromNonAlleleSpecific);
        write(jsonObject);
    }
}
