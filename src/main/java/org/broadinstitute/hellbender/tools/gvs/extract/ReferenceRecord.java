package org.broadinstitute.hellbender.tools.gvs.extract;

import htsjdk.samtools.util.Locatable;
import org.apache.avro.generic.GenericRecord;
import org.apache.avro.util.Utf8;
import org.apache.commons.lang.builder.ReflectionToStringBuilder;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.gvs.common.ChromosomeEnum;
import org.broadinstitute.hellbender.tools.gvs.common.GQStateEnum;
import org.broadinstitute.hellbender.tools.gvs.common.SchemaUtils;

import java.util.Set;

import static org.broadinstitute.hellbender.tools.gvs.common.ChromosomeEnum.*;

public class ReferenceRecord implements Locatable, Comparable<ReferenceRecord> {

    private static ChromosomeEnum chromosome;
    private final int position; // No chromosome encoded here so fits in an int for humans, dogs, cats.
    private final short length;
    private final int sampleId;
    private final short stateOrdinal;
    // Debug code -- why are there chr1's in a chr20 / chrX / chrY extract??
    private static final Set<ChromosomeEnum> expectedChromosomes = Set.of(chr20, chrX, chrY);

    public ReferenceRecord(GenericRecord genericRecord) {
        Long location = (Long) genericRecord.get(SchemaUtils.LOCATION_FIELD_NAME);
        this.position = SchemaUtils.decodePosition(location);
        setChromosome(location);
        this.sampleId = Math.toIntExact((Long) genericRecord.get(SchemaUtils.SAMPLE_ID_FIELD_NAME));

        this.length = (short) Math.toIntExact((Long) genericRecord.get("length"));
        String stringState = ((Utf8) genericRecord.get(SchemaUtils.STATE_FIELD_NAME)).toString();
        this.stateOrdinal = (short) GQStateEnum.fromValue(stringState).ordinal();
    }

    public ReferenceRecord(long location, long sampleId, int length, String stringState) {
        this.position = SchemaUtils.decodePosition(location);
        this.length = (short) length;
        setChromosome(location);
        this.sampleId = (int) sampleId;
        this.stateOrdinal = (short) GQStateEnum.fromValue(stringState).ordinal();
    }

    private static void setChromosome(Long location) {
        ChromosomeEnum thisChromosome = SchemaUtils.decodeChromosome(location);
        if (!expectedChromosomes.contains(thisChromosome)) {
            throw new GATKException("saw unexpected chromosome " + thisChromosome);
        }
        if (chromosome == null) {
            chromosome = thisChromosome;
        }
    }

    @Override
    public String getContig() { return chromosome.getContigName(); }

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
        final int positionResult = this.position - o.position;
        if (positionResult != 0) {
            return positionResult;
        }

        return this.sampleId - o.sampleId;
    }
}
