package org.broadinstitute.hellbender.tools.gvs.extract;

import htsjdk.samtools.util.Locatable;
import org.apache.avro.generic.GenericRecord;
import org.apache.commons.lang.builder.ReflectionToStringBuilder;
import org.broadinstitute.hellbender.tools.gvs.common.SchemaUtils;

public class ExtractCohortLiteFilterRecord implements Locatable {

    private final Long location;
    private final String contig;
    private final Integer start;
    private final Integer end;

    private final Double calibrationSensitivity;
    private final String yng;
    private final String refAllele;
    private final String altAllele;

    public ExtractCohortLiteFilterRecord(GenericRecord genericRecord) {
        this.location = Long.parseLong(genericRecord.get(SchemaUtils.LOCATION_FIELD_NAME).toString());
        this.contig = SchemaUtils.decodeContig(location);
        this.start = SchemaUtils.decodePosition(location);
        this.end = start;

        // TODO - this is because the table contains unfiltered records.
        Object value = genericRecord.get(SchemaUtils.CALIBRATION_SENSITIVITY);
        this.calibrationSensitivity = value != null ? Double.parseDouble(value.toString()) : 0.57;
//        this.calibrationSensitivity = Double.parseDouble(genericRecord.get("calibration_sensitivity").toString());
        this.yng = genericRecord.get(SchemaUtils.YNG_STATUS).toString();

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

    public double getCalibrationSensitivity() { return this.calibrationSensitivity; }

    public String getYng() { return this.yng; }

    public String getRefAllele() { return this.refAllele; }

    public String getAltAllele() { return this.altAllele; }

    public String toString() {
        return ReflectionToStringBuilder.toString(this);
    }
}
