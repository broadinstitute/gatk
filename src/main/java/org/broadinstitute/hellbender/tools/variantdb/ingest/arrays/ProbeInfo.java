package org.broadinstitute.hellbender.tools.variantdb.ingest.arrays;

import com.google.cloud.bigquery.FieldValue;
import com.google.cloud.bigquery.FieldValueList;
import org.broadinstitute.hellbender.utils.bigquery.QueryAPIRowReader;

import java.util.HashMap;
import java.util.Map;

public class ProbeInfo {
    public final long probeId;
    public final String name;
    public final String refBuild;
    public final String contig;
    public final String position;
    public final String ref;
    public final String alleleA;
    public final String alleleB;
    public final String flag;

    public ProbeInfo(FieldValueList fields) {
        this.probeId = fields.get(ProbeInfoSchema.PROBE_ID).getLongValue();
        this.name = fields.get(ProbeInfoSchema.NAME).getStringValue();
        FieldValue tmp = fields.get(ProbeInfoSchema.REF);
        this.ref = tmp.isNull() ? null : tmp.getStringValue();

        tmp = fields.get(ProbeInfoSchema.ALLELE_A);
        this.alleleA = tmp.isNull() ? null : tmp.getStringValue();

        tmp = fields.get(ProbeInfoSchema.ALLELE_B);
        this.alleleB = tmp.isNull() ? null : tmp.getStringValue();

        // this is specific for ingest. when we implement search or extract, we will probably need a flag
        this.refBuild = null;
        this.contig = null;
        this.position = null;
        this.flag = null;
    }

    public static Map<String, ProbeInfo> createProbeDataForIngest(QueryAPIRowReader reader) {
        Map<String, ProbeInfo> data = new HashMap<>();

        while (reader.hasNext()) {
            ProbeInfo probe = new ProbeInfo(reader.next());
            data.put(probe.name, probe);
        }
        return data;
    }

}
