package org.broadinstitute.hellbender.tools.gvs.extract;

import htsjdk.samtools.util.Locatable;
import org.apache.avro.generic.GenericRecord;
import org.apache.avro.util.Utf8;
import org.apache.commons.lang.builder.ReflectionToStringBuilder;
import org.broadinstitute.hellbender.tools.gvs.common.ChromosomeEnum;
import org.broadinstitute.hellbender.tools.gvs.common.GQStateEnum;
import org.broadinstitute.hellbender.tools.gvs.common.SchemaUtils;

public class ReferenceRecord implements Locatable, Comparable<ReferenceRecord> {

    private final ChromosomeEnum chromosome;
    private final int position; // No chromosome encoded here so fits in an int for humans, dogs, cats.
    private final short length;
    private final int sampleId;
    private final short stateOrdinal;

    public ReferenceRecord(GenericRecord genericRecord) {
        Long location = (Long) genericRecord.get(SchemaUtils.LOCATION_FIELD_NAME);
        this.position = SchemaUtils.decodePosition(location);
        this.chromosome = SchemaUtils.decodeChromosome(location);
        this.sampleId = Math.toIntExact((Long) genericRecord.get(SchemaUtils.SAMPLE_ID_FIELD_NAME));

        this.length = (short) Math.toIntExact((Long) genericRecord.get("length"));
        String stringState = ((Utf8) genericRecord.get(SchemaUtils.STATE_FIELD_NAME)).toString();
        this.stateOrdinal = (short) GQStateEnum.fromValue(stringState).ordinal();
    }

    public ReferenceRecord(long location, long sampleId, int length, String stringState) {
        this.position = SchemaUtils.decodePosition(location);
        this.length = (short) length;
        this.chromosome = SchemaUtils.decodeChromosome(location);
        this.sampleId = (int) sampleId;
        this.stateOrdinal = (short) GQStateEnum.fromValue(stringState).ordinal();
    }

    @Override
    public String getContig() { return this.chromosome.getContigName(); }

    @Override
    public int getStart() { return this.position; }

    @Override
    public int getEnd() { return this.position + this.length - 1; }

    public long getLocation() { return SchemaUtils.encodeLocation(chromosome, position); }
    public long getEndLocation() { return SchemaUtils.encodeLocation(chromosome, position + length - 1); }

    public long getSampleId() { return this.sampleId; }

    public String getState() { return GQStateEnum.fromOrdinal(stateOrdinal).getValue(); }

    public String toString() {
        return ReflectionToStringBuilder.toString(this);
    }

    @Override
    public int compareTo(ReferenceRecord o) {
        final int contigResult = this.chromosome.getIndex() - o.chromosome.getIndex();
        if (contigResult != 0) {
            return contigResult;
        }

        final int positionResult = this.position - o.position;
        if (positionResult != 0) {
            return positionResult;
        }

        return this.sampleId - o.sampleId;
    }
}
