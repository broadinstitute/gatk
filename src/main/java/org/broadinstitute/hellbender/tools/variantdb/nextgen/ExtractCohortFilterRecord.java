package org.broadinstitute.hellbender.tools.variantdb.nextgen;

import com.google.cloud.bigquery.FieldValueList;
import htsjdk.samtools.util.Locatable;
import org.apache.avro.generic.GenericRecord;
import org.broadinstitute.hellbender.tools.variantdb.SchemaUtils;

import java.util.Objects;

public class ExtractCohortFilterRecord implements Locatable {

    private final Long location;
    private final String sampleName;
    private final String contig;
    private final Integer start;
    private final Integer end;

    final Double vqslod;
    final String yng;
    private final String refAllele;
    private final String altAllele;



    public ExtractCohortFilterRecord(GenericRecord genericRecord) {
        this.location = Long.parseLong(genericRecord.get(SchemaUtils.LOCATION_FIELD_NAME).toString());
        this.sampleName = genericRecord.get(SchemaUtils.SAMPLE_NAME_FIELD_NAME).toString();
        this.contig = SchemaUtils.decodeContig(location);
        this.start = SchemaUtils.decodePosition(location);
        this.end = start;

        this.vqslod = Double.parseDouble(genericRecord.get("vqslod").toString());
        this.yng = genericRecord.get("yng_status").toString();

        // TODO - are these in fact nullable in this table?
        this.refAllele = Objects.toString(genericRecord.get(SchemaUtils.REF_ALLELE_FIELD_NAME), null);
        this.altAllele = Objects.toString(genericRecord.get(SchemaUtils.ALT_ALLELE_FIELD_NAME), null);
    }
//    public ExtractCohortFilterRecord(FieldValueList genericRecord) {
//        this.location = Long.parseLong(genericRecord.get(SchemaUtils.LOCATION_FIELD_NAME).toString());
//        this.sampleName = genericRecord.get(SchemaUtils.SAMPLE_NAME_FIELD_NAME).toString();
//        this.contig = SchemaUtils.decodeContig(location);
//        this.start = SchemaUtils.decodePosition(location);
//        this.end = start;
//        this.state = genericRecord.get(SchemaUtils.STATE_FIELD_NAME).toString();
//        this.refAllele = genericRecord.get(SchemaUtils.REF_ALLELE_FIELD_NAME).toString();
//        this.altAllele = genericRecord.get(SchemaUtils.ALT_ALLELE_FIELD_NAME).toString();
//
//    }
//

//    public ExtractCohortRecord(final Long inputLocation, final String inputSampleName) {
//        this.location = inputLocation;
//        this.sampleName = inputSampleName;
//        this.contig = SchemaUtils.decodeContig(location);
//        this.start = SchemaUtils.decodePosition(location);
//        this.end = start;
//    }



    @Override
    public String getContig() { return this.contig; }

    @Override
    public int getStart() { return this.start; }

    @Override
    public int getEnd() { return this.end; }

    public Long getLocation() { return this.location; }

    public String getSampleName() { return this.sampleName; }

    public Double getVqslod() { return this.vqslod; }

    public String getYng() { return this.yng; }

    public String getRefAllele() {
        return this.refAllele;
    }

    public String getAltAllele() {
        return this.altAllele;
    }

}
