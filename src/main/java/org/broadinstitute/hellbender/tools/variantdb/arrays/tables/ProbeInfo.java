package org.broadinstitute.hellbender.tools.variantdb.arrays.tables;

import com.google.cloud.bigquery.FieldValue;
import com.google.cloud.bigquery.FieldValueList;

import com.google.cloud.bigquery.TableResult;
import org.apache.avro.generic.GenericRecord;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.variantdb.SchemaUtils;
import org.broadinstitute.hellbender.utils.bigquery.StorageAPIAvroReader;
import org.broadinstitute.hellbender.utils.bigquery.TableReference;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

public class ProbeInfo {
    public final long probeId;
    public final String name;
    public final String refBuild;
    public final String contig;
    public final long position;
    public final String ref;
    public final String alleleA;
    public final String alleleB;
    public final String flag;

    public ProbeInfo(final long probeId, final String name, final String refBuild, final String contig,
            final long position, final String ref, final String alleleA, final String alleleB, final String flag) {
        this.probeId = probeId;
        this.name = name;
        this.contig = contig;
        this.position = position;
        this.ref = ref;
        this.alleleA = alleleA;
        this.alleleB = alleleB;
        this.refBuild = refBuild;
        this.flag = flag;
    }

    public ProbeInfo(final FieldValueList fields) {
        this.probeId = fields.get(ProbeInfoSchema.PROBE_ID).getLongValue();
        this.name = fields.get(ProbeInfoSchema.NAME).getStringValue();
        FieldValue tmp = fields.get(ProbeInfoSchema.REF);
        this.ref = tmp.isNull() ? null : tmp.getStringValue();

        tmp = fields.get(ProbeInfoSchema.ALLELE_A);
        this.alleleA = tmp.isNull() ? null : tmp.getStringValue();

        tmp = fields.get(ProbeInfoSchema.ALLELE_B);
        this.alleleB = tmp.isNull() ? null : tmp.getStringValue();

        // this is specific for ingest. when we implement search or extract, we will
        // probably need a flag
        this.refBuild = null;
        this.contig = null;
        this.position = -1;
        this.flag = null;
    }

    public static Map<String, ProbeInfo> createProbeDataForIngest(final TableResult tr) {
        final Iterator<FieldValueList> reader = tr.iterateAll().iterator();

        final Map<String, ProbeInfo> data = new HashMap<>();

        while (reader.hasNext()) {
            final ProbeInfo probe = new ProbeInfo(reader.next());
            data.put(probe.name, probe);
        }
        return data;
    }

    public static Map<Long, ProbeInfo> getProbeIdMapFromExport(final String probeCsvExportFile) {
        final Map<Long, ProbeInfo> probeIdMap = new HashMap<>();
        String line = "";

        // TODO: decide if we also want to support reading from a GCS bucket
        try (BufferedReader br = new BufferedReader(new FileReader(probeCsvExportFile))) {
            /// skip the header
            br.readLine();

            while ((line = br.readLine()) != null) {

                // use comma as separator
                final String[] fields = line.split(",");
                
                // ProbeId,Name,GenomeBuild,Chr,Position,Ref,AlleleA,AlleleB,build37Flag
                // 6,ilmnseq_rs9651229_F2BT,37,1,567667,,,,PROBE_SEQUENCE_MISMATCH
                final ProbeInfo p = new ProbeInfo(
                        Long.parseLong(fields[0]), // id
                        fields[1], // name
                        fields[2], // refBuild
                        fields[3], // contig
                        Long.parseLong(fields[4]), // position
                        fields[5], // ref
                        fields[6], // alleleA
                        fields[7], // alleleB
                        (fields.length >=9)?fields[8]:null  // build37Flag (optional)
                        );

                probeIdMap.put(p.probeId, p);
            }

            return probeIdMap;
        } catch (final Exception e) {
            throw new GATKException("Error processing probe CSV file", e);
        }
    }

    public static Map<Long, ProbeInfo> getProbeIdMapWithStorageAPI(String fqProbeTableName, boolean printDebugInformation, String readProjectId) {
        Map<Long, ProbeInfo> results = new HashMap<>();

        TableReference tableRef = new TableReference(fqProbeTableName, ProbeInfoSchema.PROBE_INFO_FIELDS);

        System.out.println("Beginning probe retrieval...");
        long start = System.currentTimeMillis();

        try (final StorageAPIAvroReader reader = new StorageAPIAvroReader(tableRef, readProjectId)) {
            for ( final GenericRecord row : reader ) {                
                ProbeInfo p = new ProbeInfo(
                    (Long) row.get(ProbeInfoSchema.PROBE_ID),
                    getOptionalString(row, ProbeInfoSchema.NAME),
                    getOptionalString(row, ProbeInfoSchema.REF_BUILD),
                    getOptionalString(row, ProbeInfoSchema.CONTIG),
                    (Long) row.get(ProbeInfoSchema.POSITION),
                    getOptionalString(row, ProbeInfoSchema.REF),
                    getOptionalString(row, ProbeInfoSchema.ALLELE_A),
                    getOptionalString(row, ProbeInfoSchema.ALLELE_B),
                    getOptionalString(row, ProbeInfoSchema.FLAG)
                );
                results.put(p.probeId, p);
            }
        }

        long elapsed = System.currentTimeMillis() - start;
        System.out.println("Completed probe retrieval... (" + elapsed + " ms)");

        return results;
    }

    private static String getOptionalString(GenericRecord rec, String fieldName) {
        Object o = rec.get(fieldName);
        return ( o == null ? null : o.toString());
    }

    public static Map<String, ProbeInfo> getProbeNameMap(final String probeCsvExportFile) {
        Map<String, ProbeInfo> results = new HashMap<>();
        for (final ProbeInfo pi : getProbeIdMapFromExport(probeCsvExportFile).values()) {
            results.put(pi.name, pi);
        }
        return results;
    }

    @Override
    public String toString() {
        return "ProbeInfo [contig=" + contig 
                + ", name=" + name + ", position=" + position + ", probeId=" + probeId + ", ref=" + ref + ", refBuild="
                + refBuild + ", alleleA=" + alleleA + ", alleleB=" + alleleB + ", flag=" + flag + "]";
    }


    

}
