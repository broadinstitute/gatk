package org.broadinstitute.hellbender.tools.gvs.extract;

import htsjdk.samtools.util.Locatable;
import org.apache.avro.generic.GenericRecord;
import org.apache.avro.util.Utf8;
import org.apache.commons.lang.builder.ReflectionToStringBuilder;
import org.apache.commons.lang3.StringUtils;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.gvs.common.ChromosomeEnum;
import org.broadinstitute.hellbender.tools.gvs.common.GQStateEnum;
import org.broadinstitute.hellbender.tools.gvs.common.SchemaUtils;

import javax.annotation.Nonnull;

public class ReferenceRecord implements Locatable, Comparable<ReferenceRecord> {

    private final int location; // no chr or contig so fits in an int, at least for humans, dogs, cats
    private final int endLocation; // same here as `location`
    private final int sampleId;
    private final ChromosomeEnum contig;
    private final GQStateEnum state;

    public ReferenceRecord(GenericRecord genericRecord) {
        Long location = (Long) genericRecord.get(SchemaUtils.LOCATION_FIELD_NAME);
        Pair<ChromosomeEnum, Integer> contigAndPosition = ChromosomeEnum.getContigAndPositionFromLocation(location);
        this.contig = contigAndPosition.getLeft();
        this.location = contigAndPosition.getRight();
        this.sampleId = (int) genericRecord.get(SchemaUtils.SAMPLE_ID_FIELD_NAME);

        int length = Math.toIntExact((Long) genericRecord.get("length"));
        this.endLocation = this.location + length - 1;
        String stringState = ((Utf8) genericRecord.get(SchemaUtils.STATE_FIELD_NAME)).toString();
        this.state = GQStateEnum.fromValue(stringState);
    }

    public ReferenceRecord(long location, long sampleId, int length, String stringState) {
        Pair<ChromosomeEnum, Integer> contigAndPosition = ChromosomeEnum.getContigAndPositionFromLocation(location);
        this.contig = contigAndPosition.getLeft();
        this.location = contigAndPosition.getRight();
        this.sampleId = (int) sampleId;
        this.endLocation = this.location + length - 1;

        this.state = GQStateEnum.fromValue(stringState);
    }

    @Override
    public String getContig() { return this.contig.getContigName(); }

    @Override
    public int getStart() { return this.location; }

    @Override
    public int getEnd() { return this.endLocation; }

    public long getLocation() { return this.location; }
    public long getEndLocation() { return this.endLocation; }

    public long getSampleId() { return this.sampleId; }

    public String getState() { return this.state.getValue(); }

    public String toString() {
        return ReflectionToStringBuilder.toString(this);
    }

    @Override
    public int compareTo(ReferenceRecord o) {
        final long firstPosition = this.location;
        final long secondPosition = o.location;

        final int result = Long.compare(firstPosition, secondPosition);
        if (result != 0) {
            return result;
        } else {
            final long firstSample = this.location;
            final long secondSample = o.location;
            return Long.compare(firstSample, secondSample);
        }
    }
}
