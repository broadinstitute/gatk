package org.broadinstitute.hellbender.tools.gvs.ingest.bq;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.tools.gvs.ingest.VetWriter;
import org.json.JSONObject;

import java.io.IOException;
import java.util.List;

public class VetBQWriter extends AbstractBQWriter implements VetWriter {

    private final boolean skipLoadingVqsrFields;
    private final boolean forceLoadingFromNonAlleleSpecific;

    public VetBQWriter(String projectId, String datasetName, String tableName, boolean skipLoadingVqsrFields, boolean forceLoadingFromNonAlleleSpecific) throws IOException {
        super(projectId, datasetName, tableName);
        this.skipLoadingVqsrFields = skipLoadingVqsrFields;
        this.forceLoadingFromNonAlleleSpecific = forceLoadingFromNonAlleleSpecific;
    }

    @Override
    public void write(long location, VariantContext variant, long sampleId, List<String> vrsAlleleIds) throws IOException {
        JSONObject jsonObject = createJson(location, variant, sampleId, skipLoadingVqsrFields, forceLoadingFromNonAlleleSpecific, vrsAlleleIds);
        write(jsonObject);
    }
}
