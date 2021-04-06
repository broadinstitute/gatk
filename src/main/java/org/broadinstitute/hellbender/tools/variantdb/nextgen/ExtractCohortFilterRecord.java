package org.broadinstitute.hellbender.tools.variantdb.nextgen;

import htsjdk.samtools.util.Locatable;
import org.apache.avro.generic.GenericRecord;
import org.broadinstitute.hellbender.tools.variantdb.SchemaUtils;

public class ExtractCohortFilterRecord implements Locatable {

    private final Long location;
    private final String contig;
    private final Integer start;
    private final Integer end;

    private final Double vqslod;
    private final String yng;
    private final String refAllele;
    private final String altAllele;

    public ExtractCohortFilterRecord(GenericRecord genericRecord) {
        this.location = Long.parseLong(genericRecord.get(SchemaUtils.LOCATION_FIELD_NAME).toString());
        this.contig = SchemaUtils.decodeContig(location);
        this.start = SchemaUtils.decodePosition(location);
        this.end = start;

        this.vqslod = Double.parseDouble(genericRecord.get("vqslod").toString());
        this.yng = genericRecord.get("yng_status").toString();

        this.refAllele = genericRecord.get(SchemaUtils.REF_ALLELE_FIELD_NAME).toString();
        this.altAllele = genericRecord.get(SchemaUtils.ALT_ALLELE_FIELD_NAME).toString();
    }

    @Override
    public String getContig() { return this.contig; }

    @Override
    public int getStart() { return this.start; }

    @Override
    public int getEnd() { return this.end; }

    public long getLocation() { return this.location; }

    public double getVqslod() { return this.vqslod; }

    public String getYng() { return this.yng; }

    public String getRefAllele() { return this.refAllele; }

    public String getAltAllele() { return this.altAllele; }

}
