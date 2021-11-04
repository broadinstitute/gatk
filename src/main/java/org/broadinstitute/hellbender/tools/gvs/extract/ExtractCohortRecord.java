package org.broadinstitute.hellbender.tools.gvs.extract;

import htsjdk.samtools.util.Locatable;
import org.apache.avro.generic.GenericRecord;
import org.apache.commons.lang.builder.ReflectionToStringBuilder;
import org.broadinstitute.hellbender.tools.gvs.common.SchemaUtils;

import java.util.Objects;

public class ExtractCohortRecord implements Locatable {

    private final Long location;
    private final Long sampleId;
    private final String contig;
    private final Integer start;
    private final Integer end;

    private final String state;
    private final String refAllele;
    private final String altAllele;
    private final String callGT;
    private final String callGQ;
    private final String callRGQ;
    private final String qualapprox;
    private final String asQualapprox;
    private final String callPL;

    // COHORT_FIELDS
//    public static final ImmutableSet<String> COHORT_FIELDS = ImmutableSet.of(
//            SchemaUtils.LOCATION_FIELD_NAME,
//            SchemaUtils.SAMPLE_ID_FIELD_NAME,
//            SchemaUtils.STATE_FIELD_NAME,
//            SchemaUtils.REF_ALLELE_FIELD_NAME,
//            SchemaUtils.ALT_ALLELE_FIELD_NAME,
//            SchemaUtils.CALL_GT,
//            SchemaUtils.CALL_GQ,
//            SchemaUtils.CALL_RGQ,
//            SchemaUtils.QUALapprox,
//            SchemaUtils.AS_QUALapprox,
//            SchemaUtils.CALL_PL);//, AS_VarDP);

    public ExtractCohortRecord(GenericRecord genericRecord) {
        this.location = Long.parseLong(genericRecord.get(SchemaUtils.LOCATION_FIELD_NAME).toString());
        this.sampleId = Long.parseLong(genericRecord.get(SchemaUtils.SAMPLE_ID_FIELD_NAME).toString());
        this.contig = SchemaUtils.decodeContig(location);
        this.start = SchemaUtils.decodePosition(location);
        this.end = start;

        // if this record is being constructed from the VET data, we won't have a state so we default it to 'v'
        this.state = Objects.toString(genericRecord.get(SchemaUtils.STATE_FIELD_NAME), "v");

        // the rest are nullable
        this.refAllele = Objects.toString(genericRecord.get(SchemaUtils.REF_ALLELE_FIELD_NAME), null);
        this.altAllele = Objects.toString(genericRecord.get(SchemaUtils.ALT_ALLELE_FIELD_NAME), null);
        this.callGT = Objects.toString(genericRecord.get(SchemaUtils.CALL_GT), null);
        this.callGQ = Objects.toString(genericRecord.get(SchemaUtils.CALL_GQ), null);
        this.qualapprox = Objects.toString(genericRecord.get(SchemaUtils.QUALapprox), null);
        this.asQualapprox = Objects.toString(genericRecord.get(SchemaUtils.AS_QUALapprox), null);
        this.callPL = Objects.toString(genericRecord.get(SchemaUtils.CALL_PL), null);

        // to keep callRGQ final...
        String tmpRGQ = Objects.toString(genericRecord.get(SchemaUtils.CALL_RGQ), null);

        // if we don't get RGQ from the database, we can calculate it from the PLs
        if (tmpRGQ != null) {
            this.callRGQ = tmpRGQ;
        } else if (this.callPL != null) {
            this.callRGQ = this.callPL.split(",")[0];
        } else {
            this.callRGQ = null;
        }
    }

    // constructor for single base reference
    public ExtractCohortRecord(long location, long sampleId, String state) {
        this.location = location;
        this.sampleId = sampleId;
        this.contig = SchemaUtils.decodeContig(location);
        this.start = SchemaUtils.decodePosition(location);
        this.end = start;
        this.state = state;

        this.refAllele = null;
        this.altAllele = null;
        this.callGT = null;
        this.callGQ = null;
        this.callRGQ = null;
        this.qualapprox = null;
        this.asQualapprox = null;
        this.callPL = null;
    }

    @Override
    public String getContig() { return this.contig; }

    @Override
    public int getStart() { return this.start; }

    @Override
    public int getEnd() { return this.end; }

    public long getLocation() { return this.location; }

    public Long getSampleId() { return this.sampleId; }

    public String getState() { return this.state; }

    public String getRefAllele() { return this.refAllele; }

    public String getAltAllele() { return this.altAllele; }

    public String getCallGT() { return this.callGT; }

    public String getCallGQ() { return this.callGQ; }

    public String getCallRGQ() { return this.callRGQ; }

    public String getQUALApprox() { return this.qualapprox; }

    public String getAsQUALApprox() { return this.asQualapprox; }

    public String getCallPL() { return this.callPL; }

    public String toString() {
        return ReflectionToStringBuilder.toString(this);
    }
}
