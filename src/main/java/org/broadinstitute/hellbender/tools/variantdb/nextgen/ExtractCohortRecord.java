package org.broadinstitute.hellbender.tools.variantdb.nextgen;

import com.google.cloud.bigquery.FieldValueList;
import htsjdk.samtools.util.Locatable;
import org.apache.avro.generic.GenericRecord;
import org.broadinstitute.hellbender.tools.variantdb.SchemaUtils;

import java.util.Objects;

//import com.google.common.collect.ImmutableSet;
//import java.util.Arrays;
//import java.util.Set;



public class ExtractCohortRecord implements Locatable {

    private final Long location;
    private final String sampleName;
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

    // optional fields from VQSR
    private final String vqslod;

    // COHORT_FIELDS
//    public static final ImmutableSet<String> COHORT_FIELDS = ImmutableSet.of(
//            SchemaUtils.LOCATION_FIELD_NAME,
//            SchemaUtils.SAMPLE_NAME_FIELD_NAME,
//            SchemaUtils.STATE_FIELD_NAME,
//            SchemaUtils.REF_ALLELE_FIELD_NAME,
//            SchemaUtils.ALT_ALLELE_FIELD_NAME,
//            SchemaUtils.CALL_GT,
//            SchemaUtils.CALL_GQ,
//            SchemaUtils.CALL_RGQ,
//            SchemaUtils.QUALapprox,
//            SchemaUtils.AS_QUALapprox,
//            SchemaUtils.CALL_PL);//, AS_VarDP);

//    [call_RGQ, ref, alt, call_GQ, QUALapprox, location, AS_QUALapprox, state, call_GT, call_PL, sample_name]

    public ExtractCohortRecord(GenericRecord genericRecord) {
        this.location = Long.parseLong(genericRecord.get(SchemaUtils.LOCATION_FIELD_NAME).toString());
        this.sampleName = genericRecord.get(SchemaUtils.SAMPLE_NAME_FIELD_NAME).toString();
        this.contig = SchemaUtils.decodeContig(location);
        this.start = SchemaUtils.decodePosition(location);
        this.end = start;
        this.state = genericRecord.get(SchemaUtils.STATE_FIELD_NAME).toString();

        // the rest are nullable
        this.refAllele = Objects.toString(genericRecord.get(SchemaUtils.REF_ALLELE_FIELD_NAME), null);
        this.altAllele = Objects.toString(genericRecord.get(SchemaUtils.ALT_ALLELE_FIELD_NAME), null);
        this.callGT = Objects.toString(genericRecord.get(SchemaUtils.CALL_GT), null);
        this.callGQ = Objects.toString(genericRecord.get(SchemaUtils.CALL_GQ), null);
        this.callRGQ = Objects.toString(genericRecord.get(SchemaUtils.CALL_RGQ), null);
        this.qualapprox = Objects.toString(genericRecord.get(SchemaUtils.QUALapprox), null);
        this.asQualapprox = Objects.toString(genericRecord.get(SchemaUtils.AS_QUALapprox), null);
        this.callPL = Objects.toString(genericRecord.get(SchemaUtils.CALL_PL), null);

        this.vqslod = Objects.toString(genericRecord.get(SchemaUtils.VQSLOD), null);
    }

    public ExtractCohortRecord(FieldValueList genericRecord) {
        this.location = Long.parseLong(genericRecord.get(SchemaUtils.LOCATION_FIELD_NAME).toString());
        this.sampleName = genericRecord.get(SchemaUtils.SAMPLE_NAME_FIELD_NAME).toString();
        this.contig = SchemaUtils.decodeContig(location);
        this.start = SchemaUtils.decodePosition(location);
        this.end = start;
        this.state = genericRecord.get(SchemaUtils.STATE_FIELD_NAME).toString();
        this.refAllele = genericRecord.get(SchemaUtils.REF_ALLELE_FIELD_NAME).toString();
        this.altAllele = genericRecord.get(SchemaUtils.ALT_ALLELE_FIELD_NAME).toString();
        this.callGT = genericRecord.get(SchemaUtils.CALL_GT).toString();
        this.callGQ = genericRecord.get(SchemaUtils.CALL_GQ).toString();
        this.callRGQ = genericRecord.get(SchemaUtils.CALL_RGQ).toString();
        this.qualapprox = genericRecord.get(SchemaUtils.QUALapprox).toString();
        this.asQualapprox = genericRecord.get(SchemaUtils.AS_QUALapprox).toString();
        this.callPL = genericRecord.get(SchemaUtils.CALL_PL).toString();

        this.vqslod = genericRecord.get(SchemaUtils.VQSLOD).toString();
    }

//    public ExtractCohortRecord(final Long inputLocation, final String inputSampleName) {
//        this.location = inputLocation;
//        this.sampleName = inputSampleName;
//        this.contig = SchemaUtils.decodeContig(location);
//        this.start = SchemaUtils.decodePosition(location);
//        this.end = start;
//    }



    @Override
    public String getContig() {
        return this.contig;
    }

    @Override
    public int getStart() {
        return this.start;
    }

    @Override
    public int getEnd() {
        return this.end;
    }

    public Long getLocation() {
        return this.location;
    }

    public String getSampleName() {
        return this.sampleName;
    }

    public String getState() {
        return this.state;
    }

    public String getRefAllele() {
        return this.refAllele;
    }

    public String getAltAllele() {
        return this.altAllele;
    }

    public String getCallGT() {
        return this.callGT;
    }

    public String getCallGQ() {
        return this.callGQ;
    }

    public String getCallRGQ() {
        return this.callRGQ;
    }

    public String getQUALApprox() {
        return this.qualapprox;
    }

    public String getAsQUALApprox() {
        return this.asQualapprox;
    }

    public String getCallPL() {
        return this.callPL;
    }

    public String getVqslod() {
        return this.vqslod;
    }
}
