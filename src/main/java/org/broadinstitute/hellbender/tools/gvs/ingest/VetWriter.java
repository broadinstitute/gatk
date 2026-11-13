package org.broadinstitute.hellbender.tools.gvs.ingest;

import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.lang3.StringUtils;
import org.json.JSONArray;
import org.json.JSONObject;

import java.io.Closeable;
import java.io.IOException;
import java.util.List;

public interface VetWriter extends Closeable {

    void write(long location, VariantContext variant, long sampleId, List<String> vrsAlleleIds) throws IOException;

    default void commitData() {
        // no-op
    }

    // Similar to create row but with types for JSON object
    default JSONObject createJson(final long location, final VariantContext variant, final long sampleId, boolean skipLoadingVqsrFields, boolean forceLoadingFromNonAlleleSpecific, final List<String> vrsAlleleIds) {
        JSONObject jsonObject = new JSONObject();
        for (final VetFieldEnum fieldEnum : VetFieldEnum.values()) {
            if (fieldEnum.equals(VetFieldEnum.location)) {
                jsonObject.put(VetFieldEnum.location.toString(), location);
            } else if (fieldEnum.equals(VetFieldEnum.sample_id)) {
                jsonObject.put(VetFieldEnum.sample_id.toString(), sampleId);
            } else if (fieldEnum.equals(VetFieldEnum.call_GQ)) {
                jsonObject.put(fieldEnum.toString(), Integer.valueOf(fieldEnum.getColumnValue(variant, forceLoadingFromNonAlleleSpecific)));
            } else {
                final String strVal = !(skipLoadingVqsrFields && fieldEnum.isVqsrSpecificField()) ? fieldEnum.getColumnValue(variant, forceLoadingFromNonAlleleSpecific) : "";
                jsonObject.put(fieldEnum.toString(), StringUtils.isEmpty(strVal) ? null : strVal);
            }
        }
        if (vrsAlleleIds != null && !vrsAlleleIds.isEmpty()) {
            jsonObject.put("vrs_allele_ids", new JSONArray(vrsAlleleIds));
        }
        return jsonObject;
    }
}
